// From a stackoverflow post
// http://stackoverflow.com/questions/22741859/deceive-the-jvm-about-the-number-of-available-cores-on-linux
//
// Setup:
// gcc -O3 -fPIC -shared -Wl,-soname,libnumcpus.so -o libnumcpus.so numcpus.c
//
// LD_PRELOAD=libnumcpus.so LIMIT_JVM_NUM_CPUS

#include <stdlib.h>
#include <unistd.h>

const char* ENV_VAR_NAME = "LIMIT_JVM_NUM_CPUS";

int JVM_ActiveProcessorCount (void) 
{
	int n_procs;
	char* cfg_val = getenv(ENV_VAR_NAME);
	if(cfg_val == NULL) {
		n_procs = sysconf(_SC_NPROCESSORS_ONLN);
	} else {
		n_procs = atoi(cfg_val);
	}
	return n_procs;
}
