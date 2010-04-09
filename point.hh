#ifndef BINACC_POINT_HH
#define BINACC_POINT_HH

// template class for holding a point

template<class T> class point
{
public:
  point() : _x(0), _y(0) {}
  point(T x, T y) : _x(x), _y(y) {}

  T x() const { return _x; }
  T y() const { return _y; }
  T& x() { return _x; }
  T& y() { return _y; }

  point& operator += (const point& p)
  { _x += p.x(); _y += p.y(); return *this; }

  point& operator *= (const T d)
  { _x *= d; _y *= d; return *this; }

  point& operator /= (const T d)
  { _x /= d; _y /= d; return *this; }

  bool operator == (const point& other) const
  {
    return _x == other._x && _y == other._y;
  }

  bool operator != (const point& other) const
  {
    return _x != other._x || _y != other._y;
  }

  // scaling operators
  point operator * ( T d ) const  { return point(_x*d, _y*d); }
  point operator / ( T d ) const  { return point(_x/d, _y/d); }

  // other operators
  point operator + (point other) const
  { return point(_x+other._x, _y+other._y); }
  point operator - (point other) const
  { return point(_x-other._x, _y-other._y); }

  // calculate distance squared from p
  T dist_sqd(const point& p) const
  {
    const T dx=x()-p.x();
    const T dy=y()-p.y();
    return dx*dx+dy*dy;
  }

private:
  T _x;
  T _y;
};

typedef point<double> point_dbl;
typedef point<int> point_int;
typedef point<unsigned short> point_ushort;

// class for comparing x coordinates of points, for sorting
template<class T> class Compare_point_x
{
public:
  int operator() (const point<T>& p1, const point<T>& p2) const
  { return p1.x() < p2.x(); }
};

// compare points in X, then Y order
template<class T> bool Compare_point_x_y (const point<T>& p1,
					  const point<T>& p2)
{
  return p1.x() < p2.x() || ( p1.x() == p2.x() && p1.y() < p2.y() );
}

// compare points in Y, then X order
template<class T> bool Compare_point_y_x (const point<T>& p1,
					  const point<T>& p2)
{
  return p1.y() < p2.y() || ( p1.y() == p2.y() && p1.x() < p2.x() );
}

#endif
