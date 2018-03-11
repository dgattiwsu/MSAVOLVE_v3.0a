#include <iostream>

using namespace std;
// definition of matrix class

template <class T>
class Mat {
private:
  int nn;
  int mm;
  T **v;
public:
  Mat();
  Mat(int n, int m);//Zero-based
  Mat(const T &a, int n, int m);//initialize to constant
  Mat(const T *a, int n, int m);//initialize to array
  Mat(const Mat &rhs);//copy constructor
  Mat & operator=(const Mat &rhs);//assignment
  Mat & operator=(const T& a);//assign a to every element
  inline T* operator[](const int i);//subscripting: pointer to row i
  inline const T* operator[](const int i) const;
  inline int nrows() const;
  inline int ncols() const;
  void print_mat();
  void print_mat(int max_col);//prints up the max columns max_col
   ~Mat();
};




template <class T>
Mat<T>::Mat():nn(0),mm(0),v(0) {}

template <class T>
Mat<T>::Mat(int n, int m): nn(n), mm(m), v(new T*[n]){
    v[0]=new T[m*n];
    //assigns pointer to the first element 
    //of an array of length m*n to v[0]
    for(int i=1;i<n;i++)
      v[i]=v[i-1]+m; //moving the pointer !
}

template <class T>
Mat<T>::Mat(const T&a, int n, int m):  nn(n), mm(m), v(new T*[n]){
  int i,j;
  v[0]=new T[m*n];
  for(i=1;i<n;i++)
    v[i]=v[i-1]+m; 
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      v[i][j]=a;
}

template <class T>
Mat<T>::Mat(const T*a, int n, int m): nn(n), mm(m), v(new T*[n]){
  int i,j;
  v[0]=new T[m*n];
  for(i=1;i<n;i++)
    v[i]=v[i-1]+m; 
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      v[i][j]=*a++;
}

template <class T>
Mat<T>::Mat(const Mat &rhs): nn(rhs.nn), mm(rhs.mm), v(new T*[nn]){
  int i,j;
  v[0]=new T[mm*nn];
  for(i=1;i<nn;i++)
    v[i]=v[i-1]+mm; 
  for(i=0;i<nn;i++)
    for(j=0;j<mm;j++)
      v[i][j]=rhs[i][j];
}

template <class T>
Mat<T> & Mat<T>::operator=(const Mat<T> &rhs){
  if(this!=&rhs){
    int i,j;
    if (nn!=rhs.nn || mm!=rhs.mm){
      if (v!=0){
	delete[] (v[0]);
	delete[] (v);
      }
      nn=rhs.nn;
      mm=rhs.mm;
      v=new T*[nn];
      v[0]=new T[mm*nn];
    }
    for(i=1;i<nn;i++)
      v[i]=v[i-1]+mm; 
    for(i=0;i<nn;i++)
      for(j=0;j<mm;j++)
	v[i][j]=rhs[i][j];
  }
  return *this;
}

template <class T>
Mat<T> & Mat<T>::operator=(const T &a){
  for(int i=0;i<nn;i++)
    for(int j=0;j<mm;j++)
      v[i][j]=a;
  return *this;
}

template <class T>
inline T* Mat<T>::operator[](const int i){
  return v[i];//defined since its a standard c array
}


template <class T>
inline const T* Mat<T>::operator[](const int i) const{
  return v[i];
}

template <class T>
inline int Mat<T>::nrows() const{
  return nn;
}

template <class T>
inline int Mat<T>::ncols() const{
  return mm;
}

template <class T>
Mat<T>::~Mat(){
  if (v!=0){
    delete[] (v[0]);
    delete[] (v);
  }
}

template <class T>
void Mat<T>::print_mat(){
  for(int i=0;i<nn;i++){
    for(int j=0;j<mm;j++){
      cout<<v[i][j]<<" ";
    }
    cout<<endl;
  }
}

template <class T>
void Mat<T>::print_mat(int max_col){
  if((max_col>nn) || (max_col>mm)){
    cerr<<"error in Mat<T>::print_mat: index out of bounds\n";
    exit(1);
  }
  for(int i=0;i<max_col;i++){
    for(int j=0;j<max_col;j++){
      cerr<<v[i][j]<<" ";
    }
    cerr<<endl;
  }
}
