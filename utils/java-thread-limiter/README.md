Java Thread Limiter
==================

A wrapper around Java applications that don't provide threadpool or process pool size control (e.g. Chemaxon).

This wrapper replaces the implementation of the JVM_ActiveProcessorCount function to return the value
of the `LIMIT_JVM_NUM_CPUS` environmental variable instead of the real processor count. This used to control certain
programs that hog resources without any override in a shared grid computing environment.

**Note**: The wrapper needs to be run at least one by someone with write access to the home directory to compile the wrappers

Example Usage(s)
-------------
     # Force molconvert to spawn only two threads
     mock-num-cpus 2 molconvert -s "c1ccccc1" mol2 -3
     
     # Use an environmental variable instead of passing in a value
     export LIMIT_JVM_NUM_CPUS=1
     mock-num-cpus molconvert -s "c1ccccc1" mol2 -3
     

System Threads and Background Processes
--------------------------------------
The wrapper also limits the number of proccesses allowed to do various system tasks to `$LIMIT_JVM_NUM_CPUS` 
by setting the following JAVA_TOOL_OPTIONS values
 * XX:CICompilerCount
 * XX:ParallelCMSThreads
 * XX:ParallelGCThreads
 
Eventually arguments will be added to allow more fine-grained control of exactly what is limited.
