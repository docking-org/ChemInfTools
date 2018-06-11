  # build an executable named myprog from myprog.c
all: best_first/best_first.c best_first/best_first.h best_first/main.c best_first_clustering/best_first.c best_first_clustering/best_first.h best_first_clustering/main.c best_first_clustering_uint16/best_first.c best_first_clustering_uint16/best_first.h best_first_clustering_uint16/main.c cal_Tc_matrix/cal_Tc_matrix.c cal_Tc_matrix/cal_Tc_matrix.h cal_Tc_matrix/main.c cal_Tc_matrix_uint16/cal_Tc_matrix_uint16.c cal_Tc_matrix_uint16/cal_Tc_matrix_uint16.h cal_Tc_matrix_uint16/main.c common/fast_tanimoto.c common/fast_tanimoto.h convert_fp_2_fp_in_16unit/main.c 
	gcc -c common/fast_tanimoto.c -o common/fast_tanimoto.o
	gcc -c best_first/best_first.c -o best_first/best_first.o
	gcc -c best_first/main.c -std=c99 -o best_first/main.o
	gcc -o best_first/best_first_clustering best_first/main.o common/fast_tanimoto.o best_first/best_first.o -lm
	gcc -c best_first_clustering/best_first.c -o best_first_clustering/best_first.o
	gcc -c best_first_clustering/main.c -std=c99 -o best_first_clustering/main.o
	gcc -o best_first_clustering/best_first_clustering best_first_clustering/main.o common/fast_tanimoto.o best_first_clustering/best_first.o -lm
	gcc -c best_first_clustering_uint16/best_first.c -o best_first_clustering_uint16/best_first.o
	gcc -c best_first_clustering_uint16/main.c -std=c99 -o best_first_clustering_uint16/main.o
	gcc -o best_first_clustering_uint16/best_first_clustering best_first_clustering_uint16/main.o common/fast_tanimoto.o best_first_clustering_uint16/best_first.o -lm
	gcc -c cal_Tc_matrix/cal_Tc_matrix.c -o cal_Tc_matrix/cal_Tc_matrix.o
	gcc -c cal_Tc_matrix/main.c -std=c99 -o cal_Tc_matrix/main.o
	gcc -o cal_Tc_matrix/cal_Tc_matrix cal_Tc_matrix/main.o common/fast_tanimoto.o cal_Tc_matrix/cal_Tc_matrix.o -lm
	gcc -c cal_Tc_matrix_uint16/cal_Tc_matrix_uint16.c -o cal_Tc_matrix_uint16/cal_Tc_matrix_uint16.o
	gcc -c cal_Tc_matrix_uint16/main.c -std=c99 -o cal_Tc_matrix_uint16/main.o
	gcc -o cal_Tc_matrix_uint16/cal_Tc_matrix_uint16 cal_Tc_matrix_uint16/main.o common/fast_tanimoto.o cal_Tc_matrix_uint16/cal_Tc_matrix_uint16.o -lm
	gcc -c convert_fp_2_fp_in_16unit/main.c -std=c99 -o convert_fp_2_fp_in_16unit/main.o
	gcc -o convert_fp_2_fp_in_16unit/convert_fp_2_fp_in_16unit convert_fp_2_fp_in_16unit/main.o common/fast_tanimoto.o  -lm


clean: 
	$(RM) */*.o 
	$(RM) best_first/best_first_clustering best_first_clustering/best_first_clustering best_first_clustering_uint16/best_first_clustering convert_fp_2_fp_in_16unit/convert_fp_2_fp_in_16unit cal_Tc_matrix_uint16/cal_Tc_matrix_uint16 cal_Tc_matrix/cal_Tc_matrix 


