
#ifndef ARRAYS_HPP_INCLUDED
#define ARRAYS_HPP_INCLUDED

#include <iostream>
#include <exception>

template <typename T> 
class ARRAY_3D
{
//---------------------------------
private:
    size_t siz;
    T* ptr;
    
    unsigned int d1,d2,d3;
    
//---------------------------------
public:
    ARRAY_3D(unsigned int _d1, unsigned int _d2, unsigned int _d3) : d1(_d1), d2(_d2), d3(_d3)
    {
//         T = nullptr;
        siz = d1*d2*d3;
        
        try
        {
            ptr = new T[siz];
        }
        catch(std::exception& e)
        {
            std::cerr << "Error while allocating internal memory for an ARRAY_3D : " << e.what() << std::endl;
        }
    
        for (unsigned int i = 0; i < siz; i++)
            ptr[i] = (T)0.0;
    }
    
    ~ARRAY_3D()
    {
        delete[] ptr;
    }

    T& operator() (unsigned int a, unsigned int b, unsigned int c)
    {
        unsigned int offset = c + b*d3 + a*d2*d3;
        return ptr[offset];
    }
    
    T  operator() (unsigned int a, unsigned int b, unsigned int c) const
    {
        unsigned int offset = c + b*d3 + a*d2*d3;
        return ptr[offset];
    }
    
    void dump()
    {
        for (unsigned int i = 0; i < siz; i++)
            std::cout << "Dump of value at rank "<< i << " : " << ptr[i] << std::endl;
    }
    
    T sum()
    {
        T l_sum = (T)0.0;
        for (unsigned int i = 0; i < siz; i++)
            l_sum += ptr[i];
        return l_sum;
    }
    
    void normalise()
    {
        T l_sum = this->sum();
        for (unsigned int i = 0; i < siz; i++)
           ptr[i] /= l_sum; 
    }
    
};

//---------------------------------------------------------------------------------------------------

template <typename T> 
class ARRAY_2D
{
//---------------------------------
private:
    size_t siz;
    T* ptr;
    unsigned int d1,d2;
    
//---------------------------------
public:
    ARRAY_2D(unsigned int _d1, unsigned int _d2) : d1(_d1), d2(_d2)
    {
//         T = nullptr;
        siz = d1*d2;
        
        try
        {
            ptr = new T[siz];
        }
        catch(std::exception& e)
        {
            std::cerr << "Error while allocating internal memory for an ARRAY_2D : " << e.what() << std::endl;
        }
    
        for (unsigned int i = 0; i < siz; i++)
            ptr[i] = (T)0.0;
    }
    
    ~ARRAY_2D()
    {
        delete[] ptr;
    }

    T& operator() (unsigned int a, unsigned int b)
    {
        unsigned int offset = b + a*d2;
        return ptr[offset];
    }
    
    T  operator() (unsigned int a, unsigned int b) const
    {
        unsigned int offset = b + a*d2;
        return ptr[offset];
    }
    
    void dump()
    {
        for (unsigned int i = 0; i < siz; i++)
            std::cout << "Dump of value at rank "<< i << " : " << ptr[i] << std::endl;
    }
    
    T sum()
    {
        T l_sum = (T)0.0;
        for (unsigned int i = 0; i < siz; i++)
            l_sum += ptr[i];
    }
    
    T normalise()
    {
        T l_sum = this->sum();
        for (unsigned int i = 0; i < siz; i++)
           ptr[i] /= l_sum; 
    }
    
};

/*
 * For testing :
 * 
int main(int argc, char* argv[])
{
    ARRAY_3D<double> a(2,2,2);
    
    double t=0.0;
    for(int i=0; i<2; i++)
        for(int j=0; j<2; j++)
            for(int k=0; k<2; k++)
            {
                a(i,j,k) = t;
                t += 1.0;
            }
            
    for(int i=0; i<2; i++)
        for(int j=0; j<2; j++)
            for(int k=0; k<2; k++)
                std::cout << a(i,j,k) << std::endl;
            
    a.dump();
    
    std::cout << "sum is : " << a.sum() << std::endl;
    
    a.normalise();
    a.dump();
    
    std::cout << "sum is : " << a.sum() << std::endl;
    
}
*/

#endif // ARRAYS_HPP_INCLUDED