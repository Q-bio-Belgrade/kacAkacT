clear all

%panel makes predictions of chaneges in steady state levels of KacA, KacT and kacAT mRNA upon changes of KacA degradation rate and the complex binding affinity;
%plots Figure 8B and 8C

kb = 0 ;
[ T_vec_0 , A_vec_0 , m_vec_0 ] = T_A_kb(kb)  ;

kb = 1 ;
[ T_vec_1 , A_vec_1 , m_vec_1 ] = T_A_kb(kb)  ;

kb = 10 ;
[ T_vec_10 , A_vec_10 , m_vec_10 ] = T_A_kb(kb)  ;

kb = 100 ;
[ T_vec_100 , A_vec_100 , m_vec_100 ] = T_A_kb(kb)  ;

kb = 1000 ;
[ T_vec_1000 , A_vec_1000 , m_vec_1000 ] = T_A_kb(kb)  ;

kb = 10000 ;
[ T_vec_10000 , A_vec_10000 , m_vec_10000 ] = T_A_kb(kb)  ;


subplot( 1 , 2 , 1)

plot( T_vec_0 , A_vec_0 , 'r' ,  T_vec_1 , A_vec_1 , 'k' , T_vec_10 , A_vec_10 , 'g' ,  T_vec_100 , A_vec_100 , 'b' , T_vec_1000 , A_vec_1000 , 'c' , T_vec_10000 , A_vec_10000 , 'm' )

subplot( 1 , 2 , 2)

plot( T_vec_0 , m_vec_0 , 'r' ,  T_vec_1 , m_vec_1 , 'k' , T_vec_10 , m_vec_10 , 'g' ,  T_vec_100 , m_vec_100 , 'b' , T_vec_1000 , m_vec_1000 , 'c' , T_vec_10000 , m_vec_10000 , 'm' )


function y = t_d(x,lambda,kb)
    phi=1;
    k=1;
    y = k*phi/(1+kb*(2/(1+lambda))^4*x^6) - x - 2*(2/(1+lambda))^4*x^6 ;

end


function [ T_vec , A_vec , m_vec ] = T_A_kb(kb)

    lambda_vec = 0:0.1:6  ;

    T_vec = [] ;

    A_vec = [] ;
    
    m_vec = [] ;

    for lambda = lambda_vec
    
        fun = @(x) t_d(x,lambda,kb);
        T = fzero(fun,1)    ;
        A = 2/(1+lambda)*T  ;
        phi = 1 ;
        m = phi/(1+kb*(2/(1+lambda))^4*T^6) ;
    
        T_vec = [ T_vec T ] ;
    
        A_vec = [ A_vec A ] ;
        
        m_vec = [ m_vec m ] ;
    
    end

end 

