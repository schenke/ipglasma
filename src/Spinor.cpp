#include "Spinor.h"

//constructor if just dimension is given
Spinor::Spinor(int n)
{
  ndim = n;
  nn = ndim;
  e = new complex<double> [nn];
  for(int i=0; i<nn; i++) e[i] = complex<double>(0.0,0.0);
}

Spinor::Spinor(int n, complex<double> a, complex<double> b, complex<double> c)
{
  if (n!=3)
    {
      cout << "error: this constructor [Spinor::Spinor(int n, complex<double> a, complex<double> b, complex<double> c)] only works for SU(3)... exiting." << endl;
      exit(1);
    }
  ndim = n;
  nn = ndim;
  e = new complex<double> [nn];
  e[0] = a;
  e[1] = b;
  e[2] = c;
}

Spinor::Spinor(int n, complex<double> a, complex<double> b)
{
  if (n!=2)
    {
      cout << "error: this constructor [Spinor::Spinor(int n, complex<double> a, complex<double> b)] only works for SU(2)... exiting." << endl;
      exit(1);
    }
  ndim = n;
  nn = ndim;
  e = new complex<double> [nn];
  e[0] = a;
  e[1] = b;
}

Spinor Spinor::normalize()
{
  return (*this)/(*this).norm();
}

Spinor Spinor::GramSchmidt(const Spinor& a)
{
  Spinor a1(ndim), c(ndim);
  a1 = a;
  a1 = a1.normalize();
  c = (*this) - (a1 * (*this)) * a1;
  return c.normalize();
}


//operators:
complex<double> operator * (const Spinor& a, const Spinor& b) // inner product
{
  int n = a.getNDim();
  complex<double> c=0.;
  for(int i=0; i<n; i++)
    {
      c += conj(a(i))*b(i);
    }
  return c;
}

//-
Spinor operator - (const Spinor& a, const Spinor& b)
{
  Spinor aa(a.getNDim());
  for(int i=0; i<a.getNN(); i++) aa.set(i,a(i)-b(i));
  return aa;
}

//+
Spinor operator + (const Spinor& a, const Spinor& b)
{
  Spinor aa(a.getNDim());
  for(int i=0; i<a.getNN(); i++) aa.set(i,a(i)+b(i));
  return aa;
}

//* multiply by a real scalar
Spinor operator * (const Spinor& a, const double s)
{
  Spinor aa(a.getNDim());
    for(int i=0; i<a.getNN(); i++){
      aa.set(i,a(i) * s);
    }
    return aa;
}
Spinor operator * (const double s, const Spinor& a)
{
  Spinor aa(a.getNDim());
    for(int i=0; i<a.getNN(); i++){
      aa.set(i,a(i) * s);
    }
    return aa;
}

//* multiply by a complex number
Spinor operator * (const complex<double> s,const Spinor& a)
{
  Spinor aa(a.getNDim());
    for(int i=0; i<a.getNN(); i++){
      aa.set(i,a(i) * s);
    }
    return aa;
}

Spinor operator * (const Spinor& a, const complex<double> s)
{
  Spinor aa(a.getNDim());
    for(int i=0; i<a.getNN(); i++){
      aa.set(i,a(i) * s);
    }
    return aa;
}

// / division by scalar
Spinor operator / (const Spinor& a, const double s)
{
  Spinor aa(a.getNDim());
    for(int i=0; i<a.getNN(); i++) aa.set(i,a(i)/s);
    return aa;
}

// complex conjugated cross product
Spinor operator % (const Spinor&  a, const Spinor& b)
{
  if (a.getNDim() == 3 )
    return Spinor(3, conj(a(1)*b(2)-a(2)*b(1)), conj(a(2)*b(0)-a(0)*b(2)),conj(a(0)*b(1)-a(1)*b(0)));
  else
    {
      cout << "[Spinor:operator %]: WARNING, not defined for N!=3" << endl;
      return 0;
    }
}


