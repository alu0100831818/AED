#include "equation_sys_t.hpp"

#include <iomanip>
#include <cmath>

namespace AED {

    equation_sys_t::equation_sys_t(void):
    A_(),
    b_()
    {}
    
    equation_sys_t::~equation_sys_t(void)
    {}
    
    void equation_sys_t::solve(matrix_t& x)
    {
        triangulariza();
        despeja(x);
    }
    
    bool equation_sys_t::verifica(const matrix_t& x, double eps) const
    {
        matrix_t Ax;
        
        Ax.multiplica(A_, x);
        
        bool iguales = true;
        for(int i = 1; (i <= Ax.get_m()) && iguales; i++){
            
            if (fabs(Ax(i, 1)- b_(i, 1)) > eps)
                iguales = false;
        }
        
        return iguales;
    }
    
    ostream& equation_sys_t::to_stream(ostream& os) const
    {
        for(int i = 1; i <= A_.get_m(); i ++){
            for(int j = 1; j <= A_.get_n(); j++)
                os << setw(10) << fixed << setprecision(4) << A_(i,j)<< " "; 
            os << " | " << setw(10) << fixed << setprecision(4) << b_(i,1)<< endl; 
        }
        return os;    
    }
    
    istream& equation_sys_t::from_stream(istream& is)
    {
        is >> A_ >> b_;
        
        return is;
    }
    
    void equation_sys_t::triangulariza(void)
    {
#ifdef _DEPURANDO_   
            clog << endl;
            clog << "--- TRIANGULARIZACIÓN de la matriz ---"<< endl;
            to_stream(clog);
            clog << endl;
            clog << "--- COMIENZO ---" << endl;
#endif        
        
        double aux;
        for(int k=1;k<A_.get_n();k++){
              for(int i=k+1;i<=A_.get_m();i++){
                if(fabs(A_(k,k))>0.0000){
                     aux=-A_(i,k)/A_(k,k);
                     A_(i,k)+=A_(k,k)*(aux);
                     for(int j=k+1;j<=A_.get_n();j++){
                             A_(i,j)+=A_(k,j)*(aux);
                     }
                      b_(i,1)+=b_(k,1)*(aux); 
                    
                 }
              }
               
#ifdef _DEPURANDO_            
            clog << "Triangulizando. Etapa "<< k << " -> "<< endl;
            to_stream(clog);
            clog << endl;
            
            
#endif
        }
        
#ifdef _DEPURANDO_            
            clog << "--- FIN ---" << endl;
#endif          
 }
    
    void equation_sys_t::despeja(matrix_t& x) const
    {
            x.resize_matrix(A_.get_n(),1);
            double aux=0.0;
            double aux1=1.0;
            for(int i=1;i<=A_.get_n();i++){
                aux1*=A_(i,i);
            }
            if(fabs(aux1)>0.0){
                     for(int i=A_.get_m();i>=1;i--){
                        for(int k=A_.get_n();k>=1;k--){
                             if(A_(i,k)!=A_(i,i))
                                aux+=A_(i,k)*x(k,1);
                         }
                         x(i,1)=(b_(i,1)-aux)/A_(i,i);
                         aux=0.0;
                     }
            }
   }

    void equation_sys_t::solut_inverse(matrix_t& x){
        /*
        A*x=B
        x=A^(-1)*B
        A^(-1)=(1/det(A))*(adj(A)^t)
        */
        if(A_.determinante()!=0){
            cout <<"determinante de la matriz: "<< A_.determinante()<< endl;
            
            matrix_t aux3;
            A_.obtener_menor(2,2,aux3);
            //MATRIZ INVERSA.
            double det_inv=(1/A_.determinante());
            A_.determinante_adj(x);
            matrix_t aux;
            aux.multiplica(x,b_);
            aux.multiplica2(det_inv);
            cout << "Solucion al sistema de ecuaciones (Metodo de Jordan): "<< endl;
            cout << aux;
            //VERFIFICACION DEL SISTEMA.
            if (verifica(aux,0.000))
                cout << "Solución correcta..."<< endl;
        }
        else{
            cerr << "No existe la matriz inversa debido a que el determinante de la matriz es 0."<< endl;
            cerr << "y por tanto no se puede resolver el sistema de ecuaciones."<< endl;
        }
    }
}

ostream&  operator<<(ostream& os, const AED::equation_sys_t& M) 
{
    return M.to_stream(os);
}

istream&  operator>>(istream& is, AED::equation_sys_t& M)
{
    return M.from_stream(is);
}