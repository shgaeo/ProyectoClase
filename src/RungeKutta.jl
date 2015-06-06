module RungeKutta
using Intervalos

export RungeKutta2, RungeKutta4, zetacero, BuscarTBuena

    
    #Método genérico de Runge-Kutta 2 para resolver ecuaciones de la forma x'=dx/dt=f(t,x)
    function RungeKutta2(t0,x0,h,n,f::Function)
        tf=t0+h*n
        t=linspace(t0,tf,n+1)

        x=similar(t); 
        x[1]=x0;

        for i in 2:n+1
            K1=h*f(t[i-1],x[i-1])   
            K2=h*f(t[i],x[i-1]+K1)
            x[i]=x[i-1]+0.5*(K1+K2)   #es la expresión general, depende de la función definida abajo
            #println(temp)
        end 
        
        
        return x
    end




    #Método genérico de Runge-Kutta 4 para resolver ecuaciones de la forma x'=dx/dt=f(t,x)
    function RungeKutta4(t0,x0,h,n,f::Function)

        tf=t0+h*n
        t=linspace(t0,tf,n+1)

        x=similar(t); 
        x[1]=x0;


        for i in 2:n+1
            K1=h*f(t[i-1]    ,x[i-1])   
            K2=h*f(t[i-1]+h/2,x[i-1]+K1/2)
            K3=h*f(t[i-1]+h/2,x[i-1]+K2/2)
            K4=h*f(t[i-1]+h  ,x[i-1]+K1)

            x[i]=x[i-1]+(1/6)*(K1+2*K2+2*K3+K4)   #es la expresión general, depende de la función definida abajo
        end 
        
        return x

    end



    function zetacero(x0,h,n,F) 
        #hasta ahora siempre t0=0
        t0=0.0;
        x_RK=RungeKutta4(t0,x0,h,n,F);
        zz=Intervalo(minimum(x_RK),maximum(x_RK)); 
        return(Intervalo(zz.a-0.001*(zz.b-zz.a) , zz.b+0.001*(zz.b-zz.a)))
    end


    function zetacero(zz::Intervalo) 
        return(Intervalo(zz.a-0.0001*(zz.b-zz.a) , zz.b+0.0001*(zz.b-zz.a)))
    end

        
    function BuscarTBuena(x0,h,n,F::Function)    
    #Esta función busca el t que hace que la ecuación 6.9 
    #del Tucker se cumpla para Psi=z0.
    #El z0 que usamos es uno que engloba la solución hasta tf 
    #y se obtiene con  Runge Kutta a 4 orden
        tf=n*h;
        z0=zetacero(x0,h,n,F) 
        
        i=0;
        a=false;
        while a==false
            T=tf-h*i
            z1=x0+Intervalo(0,T)*F(T,z0)
            a=(z1 in z0) #la condición a la que quieras llegar
            i += 1    
        end
        return (tf-h*i)
    end
        
        
        
        
    function BuscarTBuena(x0,t0,h,n,F::Function,z0::Intervalo)    
        tf=n*h+t0;
        i=0;
        a=false;
        while a==false
            T=tf-h*i+t0
            z1=x0+Intervalo(t0,T)*F(T,z0)
            a=(z1 in z0) #la condición a la que quieras llegar
            i += 1    
        end
        return (tf-h*i+t0)
    end



end
