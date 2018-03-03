#include "Cell.h"

Cell::Cell(int N)
{
  Nc = N;
  int Nc2m1 = Nc*Nc-1;
  g2mu2A = new double;
  g2mu2B = new double;
  TpA = new double;
  TpB = new double;
  U = new Matrix(Nc,1.);
  U2 = new Matrix(Nc,1.);
  Ux = new Matrix(Nc,1.);
  Uy = new Matrix(Nc,1.);
  Ux1 = new Matrix(Nc,1.);
  Uy1 = new Matrix(Nc,1.);
  Ux2 = new Matrix(Nc,1.);
  Uy2 = new Matrix(Nc,1.);
 }

Cell::~Cell()
{
  delete g2mu2A;
  delete g2mu2B;
  delete TpA;
  delete TpB;
  delete U;
  delete U2;
  delete Ux;
  delete Uy;
  delete Ux1;
  delete Ux2;
  delete Uy1;
  delete Uy2;
}

// void Cell::computeAdjointU()
// {
//   Matrix *UD;
//   UD = new Matrix(Nc);

//   *UD = *U;
//   UD->conjg();

//   complex<double> I(0,1); // this is i (0+i*1)

//   if(Nc==2)
//     {
//       complex<double> u11 = U->get(0,0);
//       complex<double> u12 = U->get(0,1);
//       complex<double> u21 = U->get(1,0);
//       complex<double> u22 = U->get(1,1);
//       complex<double> ud11 = UD->get(0,0);
//       complex<double> ud12 = UD->get(0,1);
//       complex<double> ud21 = UD->get(1,0);
//       complex<double> ud22 = UD->get(1,1);

//       UA->set(0,0.5*(u22*ud11 + u12*ud12 + u21*ud21 + u11*ud22));
//       UA->set(1,0.5*I*(u22*ud11 + u12*ud12 - u21*ud21 - u11*ud22));
//       UA->set(2,0.5*(u21*ud11 + u11*ud12 - u22*ud21 - u12*ud22));
//       UA->set(3,0.5*I*(-u22*ud11 + u12*ud12 - u21*ud21 + u11*ud22));
//       UA->set(4,0.5*(u22*ud11 - u12*ud12 - u21*ud21 + u11*ud22));
//       UA->set(5,0.5*I*(-u21*ud11 + u11*ud12 + u22*ud21 - u12*ud22));
//       UA->set(6,0.5*(u12*ud11 - u22*ud12 + u11*ud21 - u21*ud22));
//       UA->set(7,0.5*I*(u12*ud11 - u22*ud12 - u11*ud21 + u21*ud22));
//       UA->set(8,0.5*(u11*ud11 - u21*ud12 - u12*ud21 + u22*ud22));

//       //cout << *UA << endl << endl;
//     }
//   else if(Nc==3)
//     {
//       complex<double> u11 = U->get(0,0);
//       complex<double> u12 = U->get(0,1);
//       complex<double> u13 = U->get(0,2);
//       complex<double> u21 = U->get(1,0);
//       complex<double> u22 = U->get(1,1);
//       complex<double> u23 = U->get(1,2);
//       complex<double> u31 = U->get(2,0);
//       complex<double> u32 = U->get(2,1);
//       complex<double> u33 = U->get(2,2);
//       complex<double> ud11 = UD->get(0,0);
//       complex<double> ud12 = UD->get(0,1);
//       complex<double> ud13 = UD->get(0,2);
//       complex<double> ud21 = UD->get(1,0);
//       complex<double> ud22 = UD->get(1,1);
//       complex<double> ud23 = UD->get(1,2);
//       complex<double> ud31 = UD->get(2,0);
//       complex<double> ud32 = UD->get(2,1);
//       complex<double> ud33 = UD->get(2,2);
      
//       UA->set(0,0.5*(u22*ud11 + u12*ud12 + u21*ud21 + u11*ud22));
//       UA->set(1,0.5*I*(u22*ud11 + u12*ud12 - u21*ud21 - u11*ud22));
//       UA->set(2,0.5*(u21*ud11 + u11*ud12 - u22*ud21 - u12*ud22));
//       UA->set(3,0.5*(u23*ud11 + u13*ud12 + u21*ud31 + u11*ud32));
//       UA->set(4,0.5*I*(u23*ud11 + u13*ud12 - u21*ud31 - u11*ud32));
//       UA->set(5,0.5*(u23*ud21 + u13*ud22 + u22*ud31 + u12*ud32));
//       UA->set(6,0.5*I*(u23*ud21 + u13*ud22 - u22*ud31 - u12*ud32));
//       UA->set(7,0.2886751345948129*(u21*ud11 + u11*ud12 + u22*ud21 + u12*ud22 - 2.*u23*ud31 - 2.*u13*ud32));// 1/(2\sqrt(3)) = 0.2886751345948129
//       UA->set(8,0.5*I*(-u22*ud11 + u12*ud12 - u21*ud21 + u11*ud22));
//       UA->set(9,0.5*(u22*ud11 - u12*ud12 - u21*ud21 + u11*ud22));
//       UA->set(10,0.5*I*(-u21*ud11 + u11*ud12 + u22*ud21 - u12*ud22));
//       UA->set(11,0.5*I*(-u23*ud11 + u13*ud12 - u21*ud31 + u11*ud32));
//       UA->set(12,0.5*(u23*ud11 - u13*ud12 - u21*ud31 + u11*ud32));
//       UA->set(13,0.5*I*(-u23*ud21 + u13*ud22 - u22*ud31 + u12*ud32));
//       UA->set(14,0.5*(u23*ud21 - u13*ud22 - u22*ud31 + u12*ud32));
//       UA->set(15,0.2886751345948129*I*(-u21*ud11 + u11*ud12 - u22*ud21 + u12*ud22 + 2.*u23*ud31-2.*u13*ud32));
//       UA->set(16,0.5*(u12*ud11 - u22*ud12 + u11*ud21 - u21*ud22));
//       UA->set(17,0.5*I*(u12*ud11 - u22*ud12 - u11*ud21 + u21*ud22));
//       UA->set(18,0.5*(u11*ud11 - u21*ud12 - u12*ud21 + u22*ud22));
//       UA->set(19,0.5*(u13*ud11 - u23*ud12 + u11*ud31 - u21*ud32));
//       UA->set(20,0.5*I*(u13*ud11 - u23*ud12 - u11*ud31 + u21*ud32));
//       UA->set(21,0.5*(u13*ud21 - u23*ud22 + u12*ud31 - u22*ud32));
//       UA->set(22,0.5*I*(u13*ud21 - u23*ud22 - u12*ud31 + u22*ud32));
//       UA->set(23,0.2886751345948129*(u11*ud11 - u21*ud12 + u12*ud21 - u22*ud22 - 2.*u13*ud31 + 2.*u23*ud32));
//       UA->set(24,0.5*(u32*ud11 + u12*ud13 + u31*ud21 + u11*ud23));
//       UA->set(25,0.5*I*(u32*ud11 + u12*ud13 - u31*ud21 - u11*ud23));
//       UA->set(26,0.5*(u31*ud11 + u11*ud13 - u32*ud21 - u12*ud23));
//       UA->set(27,0.5*(u33*ud11 + u13*ud13 + u31*ud31 + u11*ud33));
//       UA->set(28,0.5*I*(u33*ud11 + u13*ud13 - u31*ud31 - u11*ud33));
//       UA->set(29,0.5*(u33*ud21 + u13*ud23 + u32*ud31 + u12*ud33));
//       UA->set(30,0.5*I*(u33*ud21 + u13*ud23 - u32*ud31 - u12*ud33));
//       UA->set(31,0.2886751345948129*(u31*ud11 + u11*ud13 + u32*ud21 + u12*ud23 - 2.*u33*ud31 - 2.*u13*ud33));
//       UA->set(32,0.5*I*(-u32*ud11 + u12*ud13 - u31*ud21 + u11*ud23));
//       UA->set(33,0.5*(u32*ud11 - u12*ud13 - u31*ud21 + u11*ud23));
//       UA->set(34,0.5*I*(-u31*ud11 + u11*ud13 + u32*ud21 - u12*ud23));
//       UA->set(35,0.5*I*(-u33*ud11 + u13*ud13 - u31*ud31 + u11*ud33));
//       UA->set(36,0.5*(u33*ud11 - u13*ud13 - u31*ud31 + u11*ud33));
//       UA->set(37,0.5*I*(-u33*ud21 + u13*ud23 - u32*ud31 + u12*ud33));
//       UA->set(38,0.5*(u33*ud21 - u13*ud23 - u32*ud31 + u12*ud33));
//       UA->set(39,0.2886751345948129*I*(-u31*ud11 + u11*ud13 - u32*ud21 + u12*ud23 + 2.*u33*ud31 - 2.*u13*ud33));
//       UA->set(40,0.5*(u32*ud12 + u22*ud13 + u31*ud22 + u21*ud23));
//       UA->set(41,0.5*I*(u32*ud12 + u22*ud13 - u31*ud22 - u21*ud23));
//       UA->set(42,0.5*(u31*ud12 + u21*ud13 - u32*ud22 - u22*ud23));
//       UA->set(43,0.5*(u33*ud12 + u23*ud13 + u31*ud32 + u21*ud33));
//       UA->set(44,0.5*I*(u33*ud12 + u23*ud13 - u31*ud32 - u21*ud33));
//       UA->set(45,0.5*(u33*ud22 + u23*ud23 + u32*ud32 + u22*ud33));
//       UA->set(46,0.5*I*(u33*ud22 + u23*ud23 - u32*ud32 - u22*ud33));
//       UA->set(47,0.2886751345948129*(u31*ud12 + u21*ud13 + u32*ud22 + u22*ud23 - 2.*u33*ud32 - 2.*u23*ud33));
//       UA->set(48,0.5*I*(-u32*ud12 + u22*ud13 - u31*ud22 + u21*ud23));
//       UA->set(49,0.5*(u32*ud12 - u22*ud13 - u31*ud22 + u21*ud23));
//       UA->set(50,0.5*I*(-u31*ud12 + u21*ud13 + u32*ud22 - u22*ud23));
//       UA->set(51,0.5*I*(-u33*ud12 + u23*ud13 - u31*ud32 + u21*ud33));
//       UA->set(52,0.5*(u33*ud12 - u23*ud13 - u31*ud32 + u21*ud33));
//       UA->set(53,0.5*I*(-u33*ud22 + u23*ud23 - u32*ud32 + u22*ud33));
//       UA->set(54,0.5*(u33*ud22 - u23*ud23 - u32*ud32 + u22*ud33));
//       UA->set(55,0.2886751345948129*I*(-u31*ud12 + u21*ud13 - u32*ud22 + u22*ud23 + 2.*u33*ud32 - 2.*u23*ud33));
//       UA->set(56,0.2886751345948129*(u12*ud11 + u22*ud12 - 2.*u32*ud13 + u11*ud21 + u21*ud22 - 2.*u31*ud23));
//       UA->set(57,0.2886751345948129*I*(u12*ud11 + u22*ud12 - 2.*u32*ud13 - u11*ud21 - u21*ud22 + 2.*u31*ud23));
//       UA->set(58,0.2886751345948129*(u11*ud11 + u21*ud12 - 2.*u31*ud13 - u12*ud21 - u22*ud22 + 2.*u32*ud23));
//       UA->set(59,0.2886751345948129*(u13*ud11 + u23*ud12 - 2.*u33*ud13 + u11*ud31 + u21*ud32 - 2.*u31*ud33));
//       UA->set(60,0.2886751345948129*I*(u13*ud11 + u23*ud12 - 2.*u33*ud13 - u11*ud31 - u21*ud32 + 2.*u31*ud33));
//       UA->set(61,0.2886751345948129*(u13*ud21 + u23*ud22 - 2.*u33*ud23 + u12*ud31 + u22*ud32 - 2.*u32*ud33));
//       UA->set(62,0.2886751345948129*I*(u13*ud21 + u23*ud22 - 2.*u33*ud23 - u12*ud31 - u22*ud32 + 2.*u32*ud33));
//       UA->set(63,1./6.*(u11*ud11 + u21*ud12 - 2.*u31*ud13 + u12*ud21 + u22*ud22 - 2.*(u32*ud23 + u13*ud31 + u23*ud32 - 2.*u33*ud33)));
//     }
//   else
//     {
//       cout << "[Cell::computeAdjointU]: Only SU(2) and SU(3) supported. Have Nc=" << Nc << ". Exiting." << endl;
//       exit(1);
//     }
//   delete UD;
// }
