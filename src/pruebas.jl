using FactCheck, Intervalos
using RungeKutta, PrimerOrden

facts("Probando funciones para intervalos anidados")do
    f(x)=-2x
    t0=0;tf=1.5;
    x0=3.0;z0=Intervalo(0,3);
    
    #función zetasKas
    @fact zetasKas(t0,tf,x0,z0,f) => Intervalo("-6",3)
    tf=0.5
    @fact zetasKas(tf,t0,x0,z0,f) => Intervalo(0,3.0)
    
    #función kaesimaZeta
    @fact kaesimaZeta(t0,tf,x0,z0,f,1) => [Intervalo(0,3)]
    @fact kaesimaZeta(t0,tf,x0,z0,f,10) => [Intervalo(0,3)]
    tf=0.4
    z0=Intervalo(0,3.5);
    @fact length(kaesimaZeta(t0,tf,x0,z0,f,13)) => 3
    @fact in(kaesimaZeta(t0,tf,x0,z0,f,10)[2],kaesimaZeta(t0,tf,x0,z0,f,10)[1]) => true
    @fact in(kaesimaZeta(t0,tf,x0,z0,f,10)[3],kaesimaZeta(t0,tf,x0,z0,f,10)[2]) => true
    
    #función uneInterv
    interv1=Intervalo("0.1")
    interv2=Intervalo("5.1",100000)
    interv3=Intervalo("3.1","10.1")
    @fact uneInterv([interv1,interv1],[interv2,interv3])[1] => Intervalo(interv1.a,interv2.b)
    @fact uneInterv([interv1,interv1],[interv2,interv3])[2] => Intervalo(interv1.a,interv3.b)
    
    #función conosInterv
    tf=2.0^-3; z0=Intervalo(0.5,3);
    @fact conosInterv(tf,2*tf,zetasKas(t0,tf,x0,z0,f),f) => Intervalo(1.5,2+2.0^-2+2.0^-3+2.0^-4)
end;


facts("Probando RungeKutta")do
    f(x)=-2x
    f(t,x)=-2x
    x0=3.0
    t0=0.0
    tf=1.5 #notese que 1.5=2.0^0+2.0^-1 (exactamente representable)
    n=1024 #notese que 1.5=2.0^0+2.0^-1 (exactamente representable)
    dt=tf/n
    ka=200
    pasosExtra=30;
    x_RK=RungeKutta4(t0,x0,dt,n,f)
    @fact length(x_RK) => n+1
    @fact in(Intervalo(maximum(x_RK),minimum(x_RK)),zetacero(x0,dt,n,f)) => true
    @fact BuscarTBuena(x0,dt,n,f) => less_than(tf)
    @fact BuscarTBuena(x0,dt,n,f) => greater_than(t0)
    @fact x_RK[end]-(3*exp(-2*tf)) => less_than(2.0^-5)
end;

facts("Probando función envolvente")do
    f(x)=-2x
    f(t,x)=-2x
    x0=3.0
    t0=0.0
    tf=1.5 #notese que 1.5=2.0^0+2.0^-1 (exactamente representable)
    n=1024 
    dt=tf/n #notese que 1024=2.0^10 por lo que dt=2.0^-10+2.0^-11 (exactamente representable)
    ka=200
    pasosExtra=30;
    x_RK=RungeKutta4(t0,x0,dt,n,f);
    z0=zetacero(x0,dt,n,f)
    T=BuscarTBuena(x0,dt,n,f)
    envi0=envolvente(t0,T,x0,z0,f,0);
    envi30=envolvente(t0,T,x0,z0,f,30);
    @fact length(envi0) => 1
    @fact length(envi30) => 31
    @fact in(envi30[end],envi0[end]) => true
    @fact in(x_RK[1],envi30[1]) => true
    @fact in(x_RK[T/dt],envi30[end]) => true #notese que T/dt es entero ya que dt es exactamente representable
end;

facts("Probando función NEulerMPasos")do
    t0=0.0
    tf=1.5
    n=1024
    dt=tf/n
    ordenka=200
    pasosdeintegracion=30;
    # Voy a probar algunas ODE cuya solución conozco
    
    f(x)=-2x  # solución x(t)=x0 exp(-2t)
    f(t,x)=-2x
    x0=3.0       
    eul=NEulerMPasos(f,x0,t0,n,dt,ordenka,pasosdeintegracion);
    for i=1:5
        n=int((2*rand()+1)*pasosdeintegracion/3) #paso de integración aleatorio (mayor que (1/3)pasosdeintegracion)
        @fact in(x0*exp(-2*eul[1][n])  ,  eul[3][n][1]) => true 
        # eul[1][n] es el tiempo en el n-ésimo paso de integración
        # eul[3][n][1] es la solución (primer valor de la subpartición del n_ésimo paso de integración)
    end
    
    f(x)=x  # solución x(t)=x0 exp(t)
    f(t,x)=x
    x0=Intervalo("1.1") #pruebo con condición inicial no representable
    eul=NEulerMPasos(f,x0,t0,n,dt,ordenka,pasosdeintegracion);
    for i=1:5
        n=int((2*rand()+1)*pasosdeintegracion/3) #paso de integración aleatorio (mayor que (1/3)pasosdeintegracion)
        @fact in(x0.a*exp(eul[1][n])  ,  eul[3][n][1]) => true 
        # eul[1][n] es el tiempo en el n-ésimo paso de integración
        # eul[3][n][1] es la solución (primer valor de la subpartición del n_ésimo paso de integración)
    end
     
    f(x)=-x^2  # solución x(t)=1/(t+1/x0)
    f(t,x)=-x^2
    x0=2.0
    eul=NEulerMPasos(f,x0,t0,n,dt,ordenka,pasosdeintegracion);
    for i=1:5
        n=int((2*rand()+1)*pasosdeintegracion/3) #paso de integración aleatorio (mayor que (1/3)pasosdeintegracion)
        @fact in(1/(eul[1][n]+1/x0)  ,  eul[3][n][1]) => true 
        # eul[1][n] es el tiempo en el n-ésimo paso de integración
        # eul[3][n][1] es la solución (primer valor de la subpartición del n_ésimo paso de integración)
    end
end;

