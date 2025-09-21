// Import
#include "consts.h"

// Input
#include "../../input.h"

// Functions
#include "setup.h"
#include "report.h"
#include "FFT_solve.h"
#include "DGC.h"
#include "PB.h"
#include "gets.h"

int main(void) {
	int and_hold; double freeE_old;

	SETUP();      stdout_redir_refresh();
	FFTW_SETUP(); stdout_redir_refresh();
	DGC_SETUP();  stdout_redir_refresh();
	
	REPORT_PARAMS();
	
	do {    
		iter += 1; timer_start(&t_it, &cut_it);

		timer_start(&t_dgc, &cut_dgc);
		///////////
		DGC_STEP(); stdout_redir_refresh();
		///////////
		timer_stop(&t_dgc, &cut_dgc);
		
		freeE_old = freeE;

		if ((iter % 1 == 0)){ // Currently running every time
			timer_start(&t_pb, &cut_pb);
			//////////
			PB_STEP(); stdout_redir_refresh();
			//////////
			timer_stop(&t_pb, &cut_pb);

			timer_start(&t_free, &cut_free);
			///////////////////
			freeE = get_freeE(); stdout_redir_refresh();
			///////////////////
			timer_stop(&t_free, &cut_free);

			timer_start(&t_and, &cut_and);
			/////////////////////////
			and_err = And_mix(WA, wB);
			/////////////////////////
			timer_stop(&t_and, &cut_and);
		}	
		else {
			ELEC_POLYM(); // Update rho_elec_polym if no PB_STEP
			freeE = get_freeE(); stdout_redir_refresh();
			and_err = And_mix(WA, wB);
		}
		freeDiff = fabs((freeE-freeE_old)/freeE);

		timer_stop(&t_it, &cut_it); stdout_redir_refresh();
		if (iter % 1  == 0) {
			printout();
			write_it();
			timer_report();
		}
		if (iter % 10 == 0) {
			write_W(); write_ph(); write_elec();
		}

	} while(iter<maxIter&&(and_err>Sm3||inCompMax>Sm2||freeDiff>Sm1));

	printout();
	
	write_it(); write_ph(); write_W(); write_elec();
	write_finish();
	


FILE *fp;
long int ijk; double sum=0;
fp = fopen("VOL.dat", "w");
for (ijk=0;ijk<NxNyNz;ijk++) sum+=PHA[0][ijk]*dz*dx*dy/lx/ly;
fprintf(fp, "PHA SUM                 : %10.5g\n", sum);
fprintf(fp, "EXPECTED sigma * Nm * v0: %10.5g\n", sigma_i[0]*Nm*v0);
fprintf(fp, "ERROR                   : %10.10E\n", (sum-sigma_i[0]*Nm*v0)/(sigma_i[0]*Nm*v0));
fclose(fp);


	DGC_CLEAN();
	FFTW_CLEAN();
	CLEAN();

	return 0;
}
