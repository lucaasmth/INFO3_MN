#include <sys/time.h>
#include <x86intrin.h>

#define TOP(x) gettimeofday (&x,NULL) ;

float tdiff_micro (struct timeval *start, struct timeval *end) ;

void init_flop_micro () ;

void init_flop_nano () ;

float calcul_flop_micro (char *message, unsigned int nb_operations_flottantes, float duree);

float calcul_flop_nano(char *message, int nb_operations_flottantes, unsigned long long int cycles) ;

float calcul_byte(char *message, int nb_mouvement, unsigned long long int cycles);