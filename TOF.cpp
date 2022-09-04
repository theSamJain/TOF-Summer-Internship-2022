// Time of Flight
# include <iostream>
# include <cmath>
# include <fstream>
# include <iomanip>

long double velz(long double KE, double m)
{
    long double tof, v;
    v = sqrt(2*(KE*1.60217*pow(10, -19))/(m*1.66054*pow(10, -27)));
    return(v);
}

// KE in eV
// q in C
// m in u
// theta in degrees
long double TOF(long double KE, float V[], float d[], double m, int q, long double theta)
{
    // V array ==> [V0, V1, V2 = V3, V4]
    // d array ==> [d0, d1, d2 = ~d1/2 (original position), d3, d4]
    
    long double t1, t2, t3, t4, tof, acc1, acc2, acc3, v;

    // 2ut+at^2-2d = 0  => t = -u +- root(u^2+a*2d) / a

    // Region 1
    v = velz(KE, m)*cos(theta);
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
    // std::cout<<t1<<"\n"<<t2<<"\n"<<t3<<"\n"<<t4<<std::endl;

    return(tof);
}

void Qtof(float V[], float d[])
{
    int q[] = {1, 2, 3, 4, 5};
    long double m, KE, theta;
    m = 40; // in u
    KE = 0; // in eV
    theta = 0;
    
    std::fstream newfile;
    newfile.open("A:/Projetcs & Internships/Summer Internship (TOF)/Codes/Data Files/QvTOF.dat", std::ios::out);
    for(int i = 0; i <= sizeof(q)/sizeof(q[0]) - 1; i++)
    {
        long double tof = TOF(KE, V, d, m, q[i], theta)*pow(10, 6); // microsec
        newfile<<q[i]<<"  "<<tof<<std::endl;
    }
    newfile.close();

    // std::cout<<"\nTOF (theta = "<<theta*180/3.1415269<<" deg) = "<<tof_g<<" micro seconds\n"<<std::endl;
}

void Mtof(float V[], float d[])
{
                // H, He, C, O, H2O, N2, Ar 
    long double m[] = {1, 2, 6, 8, 10, 28, 40}, KE, theta;
    int q = 1;
    KE = 0; // in eV
    theta = 0;
    
    std::fstream newfile;
    newfile.open("A:/Projetcs & Internships/Summer Internship (TOF)/Codes/Data Files/MvTOF.dat", std::ios::out);
    for(int i = 0; i <= sizeof(m)/sizeof(m[0]) - 1; i++)
    {
        long double tof = TOF(KE, V, d, m[i], q, theta)*pow(10, 6); // microsec
        newfile<<m[i]<<"  "<<tof<<std::endl;
    }
    newfile.close();
}

void MQtof(float V[], float d[])
{
                // H+, C+, C2+ O+, H2O+, N2+, N22+, Ar+, Ar2+, Ar3+
    long double m[] = {1, 6, 6, 8, 10, 28, 28, 40, 40, 40}, KE, theta;
    int q[] = {1, 1, 2, 1, 1, 1, 2, 1, 2, 3};
    KE = 0; // in eV
    theta = 0;
    
    std::fstream newfile;
    newfile.open("A:/Projetcs & Internships/Summer Internship (TOF)/Codes/Data Files/MQvTOF.dat", std::ios::out);
    for(int i = 0; i <= sizeof(m)/sizeof(m[0]) - 1; i++)
    {
        long double tof = TOF(KE, V, d, m[i], q[i], theta)*pow(10, 6); // microsec
        newfile<<double(m[i]/q[i])<<"  "<<tof<<std::endl;
    }
    newfile.close();
}

void Etof(float V[], float d[])
{
    long double KE[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, m, theta;
    int q = 1;
    m = 1; // in eV
    theta = 0;
    
    std::fstream newfile;
    newfile.open("A:/Projetcs & Internships/Summer Internship (TOF)/Codes/Data Files/EvTOF.dat", std::ios::out);
    for(int i = 0; i <= sizeof(KE)/sizeof(KE[0]) - 1; i++)
    {
        long double tof = TOF(KE[i], V, d, m, q, theta)*pow(10, 6); // microsec
        newfile<<KE[i]<<"  "<<tof<<std::endl;
    }
    newfile.close();
}

void Thetatof(float V[], float d[], int n)
{
    double theta[n], h; // Theta in Degrees
    h = 2*3.141592653/int (n-1);
    for(int i = 0; i < n; i++)
    {
        theta[i] = i*h;
    }
    long double m, KE, q;
    m = 1; // in u
    KE = 1; // in eV
    q = 1;
    
    std::fstream newfile;
    newfile.open("A:/Projetcs & Internships/Summer Internship (TOF)/Codes/Data Files/THETAvTOF.dat", std::ios::out);
    for(int i = 0; i <= sizeof(theta)/sizeof(theta[0]) - 1; i++)
    {
        long double tof = TOF(KE, V, d, m, q, theta[i])*pow(10, 6); // microsec
        newfile<<double(theta[i]*180/3.141592653)<<"  "<<tof<<std::endl;
    }
    newfile.close();
}

void sqMQtof(float V[], float d[])
{
                // H+, C+, C2+ O+, H2O+, N2+, N22+, Ar+, Ar2+, Ar3+
    long double m[] = {1, 6, 6, 8, 10, 28, 28, 40, 40, 40}, KE, theta;
    int q[] = {1, 1, 2, 1, 1, 1, 2, 1, 2, 3};
    KE = 0; // in eV
    theta = 0;
    
    std::fstream newfile;
    newfile.open("A:/Projetcs & Internships/Summer Internship (TOF)/Codes/Data Files/sqMQvTOF.dat", std::ios::out);
    for(int i = 0; i <= sizeof(m)/sizeof(m[0]) - 1; i++)
    {
        long double tof = TOF(KE, V, d, m[i], q[i], theta)*pow(10, 6); // microsec
        newfile<<pow(m[i]/q[i], 0.5)<<"  "<<tof<<std::endl;
    }
    newfile.close();
}

void gauss(float V[], float d[], int n, double mean, double std)
{
    long double E[n], h, tof_forward, tof_backward, pdfE, pdftf, pdftb;

    h = double(40)/(n-1);
    for(int i = 0; i < n; i++)
    {
        E[i] = double(20) + i*h;
        // std::cout<<E[i]<<" ";
    }
    long double m, q;
    m = 40; // in u
    q = 2;

    std::fstream newfile;
    newfile.open("A:/Projetcs & Internships/Summer Internship (TOF)/Codes/Data Files/gauss.dat", std::ios::out);
    for(int i = 0; i <= sizeof(E)/sizeof(E[0]) - 1; i++)
    {
        pdfE = exp(-(pow(((E[i]-mean)/(std)), 2))/(2)) / (std*sqrt(2*3.141592653589));
        tof_forward = TOF(E[i], V, d, m, q, 0)*pow(10, 6); // microsec
        tof_backward = TOF(E[i], V, d, m, q, 3.141592635389)*pow(10, 6); // microsec
        pdftf = 2*tof_forward*exp(-(pow(((pow(tof_forward, 2) - mean)/(std)), 2))/(2)) / (std*sqrt(2*3.141592653589));
        pdftb = 2*tof_backward*exp(-(pow((pow(tof_backward, 2)-mean)/(std), 2))/(2)) / (std*sqrt(2*3.141592653589));
        newfile<<std::setprecision(16)<<double(E[i])<<"  "<<pdfE<<" "<<tof_forward<<" "<<pdftf<<" "<<tof_backward<<" "<<pdftb<<std::endl;
    }
    newfile.close();

    std::fstream newfile2;
    double h2 = double(7)/(pow(10, 4)-1);
    newfile2.open("A:/Projetcs & Internships/Summer Internship (TOF)/Codes/Data Files/gaussTOF.dat", std::ios::out);
    for(int i = 0; i < pow(10, 4); i++)
    {
        pdftf = 2*(i*h2 + 1.5)*exp(-(pow(((pow((i*h2 + 1.5), 2) - mean)/(std)), 2))/(2)) / (std*sqrt(2*3.141592653589));
        newfile2<<std::setprecision(16)<<(i*h2+1.5)<<" "<<pdftf<<std::endl; //<<" "<<tof_backward<<" "<<pdftb<<std::endl;
    }
    newfile2.close();
}

int main()
{
    float V[] = {340, -420, -2100, -2400}; //Volts
    float d[] = {18*pow(10, -3), 8.8*pow(10, -3), 37.8*pow(10, -3), 99.8*pow(10, -3), 9.6*pow(10, -3)}; //meters

    // Qtof(V, d);
    // Mtof(V, d); 
    // Etof(V, d);
    // MQtof(V, d);
    // Thetatof(V, d, 100);
    // sqMQtof(V, d);
    // std::cout<<"done!";
    // gauss(V, d, 10000, 40, 4);

    long double mbyq[10] = {1, 7, 8, 10, 13, 18, 20, 28, 32, 40};
    for(auto i = 0; i != 10; i++)
    {
        double result = TOF(0, V, d, mbyq[i], 1, 0)/(25*pow(10, -12)) - 1339;
        std::cout<<std::setprecision(16)<<result<<std::endl;
    }

    // double result = TOF(0, V, d, 1, 1, 0)/(25*pow(10, -12)) - 1339;
    // std::cout<<std::setprecision(16)<<result<<std::endl;
}