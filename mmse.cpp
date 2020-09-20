/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#include "utilities.h"
#include "mmse.h"

MMSE::MMSE(string x, string y) : eq(true)
{
  string str;
  cout << "How many variables will your approximation have?" << endl;
  getline (cin, str);
  if (!VerifyInput(str, false))
  {
   cout << "Number of variables has to be of type unsigned int\n";
   exit(ERROR_BAD_INPUT);
  }
  stringstream(str) >> this-> numVar;

  if (x != "")
  {
    this->X.SetFromFile(x);
  }
  else
  {
    cout << "How many points would you like to enter? " << endl;
    getline (cin, str);
    if (!VerifyInput(str, false))
    {
     cout << "Number of variables has to be of type unsigned int\n";
     exit(ERROR_BAD_INPUT);
    }
    stringstream(str) >> this->X.m;
    this->X.n = 1;
    cout << "Regarding the X values, please enter: " << endl;
    this->X.SetMatrix();  
  }
  if (y != "")
  {
    this->Y.SetFromFile(y);
  }
  else
  {
    if(X.m != 0)
    {
      this->Y.m = this->X.m;
      this->Y.n = 1;
      cout << "Regarding the Y values, please enter: "<< endl;
      this->Y.SetMatrix();
    }
    else
    {
      cout << "You have to enter at least 1 coordinate for the SSME methods" << endl;
      exit (ERROR_BAD_INPUT);
    }
  }  
  
}

void MMSE::Solve()
{
  int i = 0;
  this->P.m = this->X.m;
  this->P.n = 2;
  P.Allocate();
  
  for (i = 0; i < int(P.m); i++)
  {
   this->P.matrix[i][0] = 1;
   this->P.matrix[(unsigned int) i][1] = this->X.matrix[(unsigned int) i][0];
  }
  
  this->P.Transpose(this->Pt);
  this->eq.numberVariables = this-> numVar; 
  this->Pt.Cross(this->eq.matrixA, this->P);
  this->Pt.Cross(eq.matrixB, this->Y);
  eq.Check();
  eq.Solve();
}

