#include<madness/chem/SCF.h>
#include<madness/chem/commandlineparser.h>
#include<madness/chem/molopt.h>
#include <madness/world/worldmem.h>

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && defined(HAVE_UNISTD_H)

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include<madness/chem/QCCalculationParametersBase.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include<madness/world/worldgop.h>


static inline int file_exists(const char *inpname) {
    struct stat buffer;
    int rc = stat(inpname, &buffer);
    return (rc == 0);
}

#endif

using namespace madness;


static double ttt, sss;

static void START_TIMER(World& world) {
    world.gop.fence();
    ttt = wall_time();
    sss = cpu_time();
}

static void END_TIMER(World& world, const char *msg) {
    ttt = wall_time() - ttt;
    sss = cpu_time() - sss;
    if (world.rank() == 0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
}

int main(int argc, char **argv) {

    initialize(argc, argv);

    { // limit lifetime of world so that finalize() can execute cleanly
        World world(SafeMPI::COMM_WORLD);
        START_TIMER(world);
        // Load info for MADNESS numerical routines
        startup(world, argc, argv, true);
        commandlineparser parser(argc, argv);
        std::cout.precision(6);
        SCF calc(world, parser);
        
       	
        if (world.rank() == 0) std::cout << calc.molecule.get_all_coords().flat() << "\n";	
        calc.set_protocol<3>(world, calc.param.protocol()[0]);
        MolecularEnergy E(world, calc);
        double energy = E.value(calc.molecule.get_all_coords().flat());
        if ((world.rank() == 0) and (calc.param.print_level() > 0))
            printf("final energy=%16.8f ", energy);
        E.output_calc_info_schema();

	auto ao = calc.ao;
        if(world.rank()==0){
	    std::cout << ao.size() << " atomics\n";
	}
	auto left = ao[0];
	auto op = calc.coulop; 
	auto rho = left*left;
	auto J = apply(*op, rho);
	const double integral = (left*left).inner(J);
        if(world.rank()==0){
            std::cout << "integral: " << integral << "\n";
	}

        // save the J potential
	madness::save(J, "jpot");       

        // read in gridpoints
	std::string filename = "gridpoints.txt";
	std::string line,filecontents;
	std::ifstream fs(filename.c_str());
        while (std::getline(fs, line)) filecontents += line + "\n";

        std::stringstream f(filecontents);
//position_stream_to_word(f, tag, '#', true, true);

        // read input lines
        while (std::getline(f,line)) {
             std::stringstream sline(line);	
             double x,y,z;
             sline >> x;
             sline >> y;
             sline >> z;
             madness::coord_3d R;
	     R[0]=x;
	     R[1]=y;
	     R[2]=z;
             double jx = J.eval(R);
             if(world.rank() == 0) std::cout << R << " : " << jx << "\n";	    
	} 
    }
    finalize();
    return 0; 
}
