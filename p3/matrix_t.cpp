#include "matrix_t.hpp"

#include <iomanip>
#include <cstdio>
#include <cmath>

using namespace std;

namespace AED {
  
    int matrix_t::pos(int i,int j) const
    {
    	if ((i<1)||(i>m_)||(j<1)||(j>n_)){
    		cerr << "ERROR: acceso a matriz incorrecto"<< endl;
    		return 0;
        }
    
    	return (i-1)*n_+j-1;
    }
    
    
    void matrix_t::build_matrix(void) 
    {
    	M_= new double [m_ * n_];
    	
    	if (M_==NULL)  
    		cerr << "ERROR: memoria insuficiente" << endl;
    }		
    
    void matrix_t::destroy_matrix(void)
    {
    	if (M_!=NULL){
    		delete [] M_;	
    		M_=NULL;	
    	}
    
    	m_=0;
    	n_=0;
    }
    
    void matrix_t::resize_matrix(int m,int n)
    {
    	destroy_matrix();
    	m_=m;
    	n_=n;
    	build_matrix();
    }
    
    matrix_t::matrix_t(int m,int n):
    M_(NULL),
    m_(m),
    n_(n)
    {
    	build_matrix();
    }	
    
    matrix_t::matrix_t(void):
    M_(NULL),
    m_(0),
    n_(0)
    {}	
    
    matrix_t::matrix_t(const matrix_t& M):
    M_(NULL),
    m_(M.get_m()),
    n_(M.get_n())
    {
        resize_matrix(m_, n_);
        
        for(int i = 1; i <= m_; i ++)
            for(int j = 1; j <= n_; j++)
                get_set_matrix_item(i,j) = M(i,j);
    }
    
    matrix_t::~matrix_t(void)
    {
        destroy_matrix();
    }
    
    istream&  matrix_t::from_stream(istream& is)
    {
        destroy_matrix();
        
        is >> m_ >> n_;
        
        resize_matrix(m_, n_);
        
    	for(int i=1; i <= m_; i++)
    		for(int j=1;j<=n_;j++)
    		    is >> get_set_matrix_item(i,j);
    		    
    	return is;
    }
    
    ostream&  matrix_t::to_stream(ostream& os) const
    {
        //os << setw(10) << m_ << endl;
        //os << setw(10) << n_ << endl;
        os << endl;
    
    	for(int i=1; i <= m_; i++){
    	    os<< "|";
    		for(int j=1;j<=n_;j++)
    		    os << setw(10) << fixed << setprecision(2) << get_matrix_item(i,j)<< " ";
    		
    		os << "|"<< endl;
    	}
    
    	os << endl;
    	
    	return os;
    }
    
    double matrix_t::get_matrix_item(int i,int j) const
    {
    	return M_[pos(i,j)];
    }
    
    double& matrix_t::get_set_matrix_item(int i,int j)
    {
    	return M_[pos(i,j)];
    }
    
    int matrix_t::get_m(void) const
    {
    	return m_;
    }
    
    int matrix_t::get_n(void) const
    {
    	return n_;
    }
    
	double matrix_t::operator()(int i, int j) const
	{
	    return get_matrix_item(i,j);   
	}
	
	double& matrix_t::operator()(int i, int j)
	{
	    return get_set_matrix_item(i, j);
	}
    
    void matrix_t::multiplica(const matrix_t& A, const matrix_t& B)
    {
        if(multiplicable(A, B)){
            resize_matrix(A.get_m(), B.get_n());
            
            for(int i = 1; i <= m_; i ++)
                for(int j = 1; j <=n_; j++){
                    
                    double aux = 0.0;
                    
                    for(int k = 1; k <= A.get_n(); k++)
                        aux += A(i,k) * B(k,j);
                        
                    get_set_matrix_item(i,j) = aux;
                }
        }
        else{
            cerr << "ERROR: Matrices no multiplicables"<< endl; 
        }
    }
    
    bool matrix_t::multiplicable(const matrix_t& A, const matrix_t& B)
    {
        if (A.get_n() == B.get_m())
            return true;
        else
            return false;
    }

    bool matrix_t::cuadrada(void)
    {
        if (m_ == n_)
            return true;
        else
            return false;
    }
     void matrix_t::adjunta(matrix_t& a,int i_i, int j_j){
        a.resize_matrix(m_-1,n_-1);
        int f=0;
        for(int i=1;i<=a.get_m();i++){
                if(i==i_i) f=1;
                int c=0;
                for(int j=1;j<=a.get_n();j++){
                    if(j==j_j) c=1;
                    a(i,j)=get_matrix_item(i+f,j+c);
                }
        }
    }
    
    void matrix_t:: obtener_menor(int i,int j, matrix_t& m){
        adjunta(m,i,j);
        cout <<"menor de orden "<< m_-1 <<" x "<< n_-1<<" " << endl;
        cout << m<< endl;
    }
    
    int matrix_t::signo(int i,int j){
        //return (i+j)%2?1:-1; ()
        if((i+j)%2==0) return 1;
        else return -1;
    }
    
    double matrix_t::determinante(void){
        if(m_==1) return get_matrix_item(1,1);
        else{
            double det=0;
            for(int i=1;i<=m_;i++){
                matrix_t aux;
                adjunta(aux,i,1);
                det+=get_matrix_item(i,1)*signo(i,1)*aux.determinante();
            }
            return det;
        }
    }
    void matrix_t::determinante_adj(matrix_t& x){
        x.resize_matrix(get_m(),get_n());
        if(m_==1) x(1,1)= get_matrix_item(1,1);
        else{
            for(int i=1;i<=m_;i++){
                for(int j=1;j<=m_;j++){
                    matrix_t aux;
                    adjunta(aux,i,j);
                    x(j,i)=signo(i,j)*aux.determinante();
                }
            }
            
        }
    }
    
    void matrix_t:: multiplica2(double y){
        for(int i=1;i<=get_m();i++)
            for(int j=1;j<=get_n();j++)
                get_set_matrix_item(i,j)*=y;
    }

}


ostream&  operator<<(ostream& os, const AED::matrix_t& M)
{
    return M.to_stream(os);
}

istream&  operator>>(istream& is, AED::matrix_t& M)
{
    return M.from_stream(is);
}