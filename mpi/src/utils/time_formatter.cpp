#include <string>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <math.h>

using namespace std;


void getFormattedTime(double seconds, char *formatted_string)
{
    double s(fabs(seconds));
    int h(s/3600);
    int min(s/60 - h*60);
    double sec(s - (h*60 + min)*60);
    std::ostringstream oss;
    oss<<std::setfill('0')<<std::setw(2)<<fabs(seconds)/seconds*h<<":"<<std::setw(2)<<min<<":";
    if (sec/10<1)
        oss<<"0";
    oss<<sec;
    strcpy(formatted_string,oss.str().c_str());
}