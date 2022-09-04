# include <iostream>
# include <cmath>
# include <stdlib.h>
# include <iomanip>
# include <fstream>
# include <vector>

long double TOF_Vz(long double vz, float V[], float d[], double m, int q)
{
    // V array ==> [V0, V1 = V2, V3, V4]
    // d array ==> [d0, d1, d2 = ~d1/2 (original position), d3, d4]
    
    long double t1, t2, t3, t4, tof, acc1, acc2, acc3, v;

    // 2ut+at^2-2d = 0  => t = -u +- root(u^2+a*2d) / a

    // Region 1
    v = vz;
    acc1 = -(q*1.60217*pow(10, -19)*(V[1]-V[0]))/(d[0]*m*1.66054*pow(10, -27));
    t1 = (-v + sqrt(pow(v, 2) + 2*d[1]*acc1))/acc1;
    
    // Region 2
    v += acc1*t1;
    acc2 = -(q*1.60217*pow(10, -19)*(V[2]-V[1]))/(d[2]*m*1.66054*pow(10, -27));
    t2 = (-v + sqrt(pow(v, 2) + 2*d[2]*acc2))/acc2;

    // Region 3
    v += acc2*t2;
    t3 = (d[3])/(v);
    
    // Region 4
    acc3 = -(q*1.60217*pow(10, -19)*(V[3]-V[2]))/(d[4]*m*1.66054*pow(10, -27));
    t4 = (-v + sqrt(pow(v, 2) + 2*d[4]*acc3))/acc3;

    tof = t1 + t2 + t3 + t4;

    return(tof);
}

long double func(long double x, double m, int q, long double tof)
{
    long double funcval;
    long double t_offset = 1339*25*pow(10, -12); 

    float V[] = {340, -420, -2100, -2400}; //Volts
    float d[] = {18*pow(10, -3), 8.8*pow(10, -3), 37.8*pow(10, -3), 99.8*pow(10, -3), 9.6*pow(10, -3)}; //meters


    funcval = (tof + t_offset) - TOF_Vz(x, V, d, m, q); 
    return(funcval);
}

long double bisection(long double a, long double b, long double epsilon, long double param[])
{   
    long double res, c;
    double m = param[0];
    int q = int(param[1]);
    long double tof = param[2];

    if(func(a, m, q, tof)*func(b, m, q, tof) > 0)
    {
        std::cout<<"Wrong Interval";
        std::cout<<tof<<" "<<m;
        exit(1);
    }
    
    if(func(a, m, q, tof)*func(b, m, q, tof) <= 0)
    {
        while((b-a) >= epsilon)
        {
            c = (a+b)/2;
            // std::cout<<c<<" ";
            if(func(c, m, q, tof) == 0)
            {
                break;
            }
            else if(func(a, m, q, tof)*func(c, m, q, tof) < 0)
            {
                b = c;
            }
            else
            {
                a = c;
            }
            res = c;
        }
    }
    return(res);
}

void KER_spec(long double a, long double b, long double eps)
{
    std::fstream tof_file;
    tof_file.open("A:/Projetcs & Internships/Summer Internship (TOF)/Codes/Data Files/Time.txt", std::ios::in);

    std::vector<double> tof1, tof2;
    while(!tof_file.eof())
    {
        double a, b;
        tof_file>>a>>b;
        tof1.push_back(a*25*pow(10, -12));
        tof2.push_back(b*25*pow(10, -12));
    }

    long double m1, q1, m2, q2, KER1, KER2, vz1, vz2, KER;

    m1 = 15.0;
    q1 = 1.0;
    m2 = 39.0;
    q2 = 1.0;

    std::fstream vz_file;
    vz_file.open("A:/Projetcs & Internships/Summer Internship (TOF)/Codes/Data Files/Vz_new.txt", std::ios::out);
    for(auto j = 0; j != tof1.size()-1; j++)
    {
        long double param1[] = {m1, int(q1), tof1[j]};
        long double param2[] = {m2, int(q2), tof2[j]};

        vz1 = bisection(a, b, eps, param1)/(2.18769*pow(10, 6)); // in au
        // KER1 = (m1*1823.391*pow(vz1, 2)/2)*(27.211386); // in eV
        
        vz2 = bisection(a, b, eps, param2)/(2.18769*pow(10, 6)); // in au
        // KER2 = m2*1823.391*pow(vz2, 2)/(2)*(27.211386);  // in eV
        
        // KER = KER1 + KER2;  // in eV
        
        // vz_file<<std::setprecision(13)<<vz1<<" "<<KER1<<" "<<vz2<<" "<<KER2<<" "<<KER<<std::endl;
        vz_file<<std::setprecision(13)<<vz1<<" "<<vz2<<std::endl;
    }

    vz_file.close();
}

int main()
{
    // double m = 15;
    // int q = 1;
    // long double tof = 9.96842e-006;
    // long double param[3] = {m, q, tof};
    // long double a = -pow(10, 7), b = pow(10, 7), res;
    // long double eps = pow(10, -13);
    // res = bisection(a, b, eps, param);
    // std::cout<<std::setprecision(15)<<res/(2.18769*pow(10, 6));

    long double a = -pow(10, 7), b = pow(10, 7), res;
    long double eps = pow(10, -12);
    KER_spec(a, b, eps);
    std::cout<<"done!";
}