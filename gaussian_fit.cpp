//#include <opencv2/opencv.hpp>
using namespace std;
std::vector<double> params(std::vector<double> x,std::vector<double> y) {

        double Y[201], sumx=0, sumy=0,sumxy=0, sumx2=0, ex=0, A, b;
        std::vector<double> func;

        for(int i=0; i<=200; i++) {
                Y[i] = log(y.at(i));
        }

        double mean, sigma;
        for(int i=0; i<=200; i++) {
                sumx = sumx + (x.at(i));
                sumx2=sumx2 + (x.at(i) * x.at(i));
                sumy=sumy + Y[i];
                sumxy=sumxy + (x.at(i)) * Y[i];
                mean = (sumxy / sumy);
                ex = (Y[i] * (x.at(i) - mean) * (x.at(i) - mean));
                sigma = sqrt(ex/sumy)+0.0215;
                A=((sumx2*sumy -sumx*sumxy)*1.0/(201*sumx2-sumx*sumx)*1.0);
        //      b=((201*sumxy-sumx*sumy)*1.0/(201*sumx2-sumx*sumx)*1.0);
   //           func.push_back(A*exp((-1*(x.at(i)-mean)*(x.at(i)-mean))/2*sigma*sigma));
                //func.push_back(A*exp(b));
        }

        func.push_back(A);
        func.push_back(mean);
        func.push_back(sigma);
        return func;


       // printf("\n\n The curve is Y= %fe^%fX",a,b);
}

void curve_fit(std::vector<double> x, std::vector<double> y, gnuplot_ctrl *h, int num) {
        //srand(201);

        std::vector<double> func = params(x,y);
        std::vector<double> y2;
        for(int i=0; i<x.size(); i++){
                if(num == 2 || num == 3) {
                        y2.push_back(-func[0]*exp(-1*pow((x.at(i)-func[1]), 2)/(2*pow(func[2], 2))));
                }
                else {
                        y2.push_back(func[0]*exp(-1*pow((x.at(i)-func[1]), 2)/(2*pow(func[2], 2))));
                }
        }

        double Y[201], X[201];
        /*
        double ymax = y2.at(0);
        for(int i=1; i<=200; i++) {
                if(y2.at(i) > ymax) {
                        ymax = y2.at(i);
                }

        }
        std::vector<double> xnew, ynew;
        for(int i=0; i<x.size(); i++) {
                if(y2.at(i) > ymax*1.0) {
                        xnew.push_back(x.at(i));
                        ynew.push_back(y2.at(i));
                }
        }
        */
        
        
        for(int i=0; i<=200; i++) {

                if(num == 1) {
                        Y[i] = 0.1*(log(y2.at(i))/10)+4.259178;
                }
                if(num == 2) {
                        Y[i] = 0.1*(log(y2.at(i))/10)+2.663588;
                        //+2.663588;
                }
                if(num == 3) {
                        Y[i] = 0.1*(log(y2.at(i))/10)+3.640968;
                        //+3.651068;
                }

                X[i] = x.at(i);
        }

       gnuplot_setstyle(h, "lines");
       gnuplot_plot_xy(h, X, Y, 201, "fit");

          cout << "Polynomial fit!" << endl;

    // Input values
    // **************************************************************
    size_t k = 2;                                    // Polynomial order
    bool fixedinter = false;                         // Fixed the intercept (coefficient A0)
    int wtype = 0;                                   // Weight: 0 = none (default), 1 = sigma, 2 = 1/sigma^2
    double fixedinterval = 0.;                       // The fixed intercept value (if applicable)
    double alphaval = 0.05;                          // Critical apha value

    double erry[] = {0.1, 0.3, 0.2, 0.4, 0.1, 0.3};       // Data points (err on y) (if applicable)

    // Definition of other variables
    // **************************************************************
    size_t n = 201;                                    // Number of data points (adjusted later)
    size_t nstar = 0;                                // equal to n (fixed intercept) or (n-1) not fixed
    double coefbeta[k+1];                            // Coefficients of the polynomial
    double serbeta[k+1];                             // Standard error on coefficients
    double tstudentval = 0.;                         // Student t value
    double SE = 0.;                                  // Standard error

    double **XTWXInv;                                // Matrix XTWX Inverse [k+1,k+1]
    double **Weights;                                // Matrix Weights [n,n]


    // Initialize values
    // **************************************************************
   // n = sizeof(x)/sizeof(double);
    nstar = n-1;
    if (fixedinter) nstar = n;

    cout << "Number of points: " << n << endl;
    cout << "Polynomial order: " << k << endl;
    if (fixedinter) {
        cout << "A0 is fixed!" << endl;
    } else {
        cout << "A0 is adjustable!" << endl;
    }

    /*if (k>nstar) {
        cout << "The polynomial order is too high. Max should be " << n << " for adjustable A0 ";
        cout << "and " << n-1 << " for fixed A0. ";
        cout << "Program stopped" << endl;
        return -1;
    }*/

    if (k==nstar) {
        cout << "The degree of freedom is equal to the number of points. ";
        cout << "The fit will be exact." << endl;
    }

    XTWXInv = Make2DArray(k+1,k+1);
    Weights = Make2DArray(n,n);

    // Build the weight matrix
    // **************************************************************
    CalculateWeights(erry, Weights, n, wtype);

    cout << "Weights" << endl;
    displayMat(Weights,n,n);

/*    if (determinant(Weights,n)==0.) {
        cout << "One or more points have 0 error. Review the errors on points or use no weighting. ";
        cout << "Program stopped" << endl;
        return -1;
    } */

    // Calculate the coefficients of the fit
    // **************************************************************
    PolyFit(X,Y,n,k,fixedinter,fixedinterval,coefbeta,Weights,XTWXInv);


    // Calculate related values
    // **************************************************************
    double RSS = CalculateRSS(X,Y,coefbeta,Weights,fixed,n,k+1);
    double TSS = CalculateTSS(X,Y,coefbeta,Weights,fixedinter,n,k+1);
    double R2 = CalculateR2COD(X,Y,coefbeta,Weights,fixedinter,n,k+1);
    double R2Adj = CalculateR2Adj(X,Y,coefbeta,Weights,fixedinter,n,k+1);

    if ((nstar-k)>0) {
        SE = sqrt(RSS/(nstar-k));
        tstudentval = fabs(CalculateTValueStudent(nstar-k, 1.-0.5*alphaval));
    }
    cout << "t-student value: " << tstudentval << endl << endl;

    // Calculate the standard errors on the coefficients
    // **************************************************************
    CalculateSERRBeta(fixedinter,SE,k,serbeta,XTWXInv);

    // Display polynomial
    // **************************************************************
    DisplayPolynomial(k);

    // Display polynomial coefficients
    // **************************************************************
    DisplayCoefs(k, nstar, tstudentval, coefbeta, serbeta);

    // Display statistics
    // **************************************************************
    DisplayStatistics(n,nstar,k,RSS,R2,R2Adj,SE);

    // Display ANOVA table
    // **************************************************************
    DisplayANOVA(nstar, k, TSS, RSS);

    // Write the prediction and confidence intervals
    // **************************************************************
    WriteCIBands("CIBands2.dat",X,coefbeta,XTWXInv,tstudentval,SE,n,k);

    // Display the covariance and correlation matrix
    // **************************************************************
    DisplayCovCorrMatrix(k, SE, fixedinter, XTWXInv);
}