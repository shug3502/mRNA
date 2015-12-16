#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char *argv[ ])
{
int int_is_anchored = 0;
int int_is_attached = 1;
int i = 0;
int j = 0;

float float_params[6] = {1.16, 0.80, 0.11, 0.42, 0.84, 0.58}	;
float float_time = 0;
float float_L[2] = {52, 37};
float float_nuc_rad = 10;
float float_rc_width = 1;
float float_xpos = 0;
float float_ypos = 0;
float float_alpha = 0;
float float_rr[4];
float float_tau=10;
float float_delx;
float float_theta;

char charArr25_inputValue[25];
printf("What is the end time? ");
scanf("%24s", charArr25_inputValue);
float float_end_time = atof(charArr25_inputValue);


while (float_time < float_end_time){
// Total Reaction rate 
float_alpha = float_params[2]*(1-int_is_attached) + float_params[3] + float_params[4];

//rr = rand(4,1);
float_tau = 1/float_alpha*log(1/rr[0]);
//delx = tau*v;	


float_time = float_time + float_tau;
printf("%f\n",float_time);
}




} // end main