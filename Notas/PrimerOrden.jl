module PrimerOrden
using Intervalos
using RungeKutta
using PyPlot



export zetasKas, kaesimaZeta, graficaCono,conosInterv,uneInterv,envolvente,graficaEnvolvente,ptosenvolventes,colores,graficarEulerM,NEulerMPasos,NEuler1Paso


function zetasKas(tI,tF,x0,z0::Intervalo,f::Function)  
    x0 + Intervalo(tI,tF)*f(z0)
end


function kaesimaZeta(tI,tF,x0,z0::Intervalo,f::Function,k::Int64)
    zetas=[z0]
    for i=1:k
        zk=zetasKas(tI,tF,x0,zetas[i],f)
        if zk==zetas[i]
            return zetas
        else
            push!(zetas,zk)
        end
    end
    return zetas
end





####### PARA GRAFICAR Y PARA SACAR LOS SUBINTERVALOS


function graficaCono(tI,tF,x0,zk::Intervalo,col::ASCIIString)
    PyPlot.fill_between([tI,tF],[x0,float64(zk.a)],[x0,float64(zk.b)],color=col)
    #PyPlot.xlim(min(0.9*tI,1.1*tI),max(0.9*tF,1.1*tF))
    #PyPlot.ylim(min(0.9*float64(zk.a),1.1*float64(zk.a)),max(0.9*float64(zk.b),1.1*float64(zk.b)))
    return nothing
end

graficaCono(tI,tF,x0,zk::Intervalo)=graficaCono(tI,tF,x0,zk,"lightblue")


function conosInterv(tI,tF,x0::Intervalo,f::Function)
    t=tF-tI
    cono1=x0.a + t*f(x0)    
    cono2=x0.b + t*f(x0)
    conote=uneInterv([cono1],[cono2])
    return conote[1]
end

function uneInterv(conoA::Array{Intervalo,1},conoB::Array{Intervalo,1})
    # esta función toma dos arreglos de intervalos y los une aunque no se intersecten
    if length(conoA)!=length(conoB)
        error("Deben tener mismo tamaño")
    end
    if conoA[1].a≥conoB[1].a && conoA[end].a≥conoB[end].a
        bajo=conoB
        alto=conoA
    else
        bajo=conoA
        alto=conoB
    end
    cono=Intervalo[]
    for i=1:length(conoA)
        push!(cono,Intervalo((bajo[i]).a,(alto[i]).b))
    end
    return cono
end



function graficaCono(tI,tF,x0::Intervalo,zk::Intervalo,lineas::Bool,col::ASCIIString)
    PyPlot.fill_between([tI,tF],Float64[x0.a,zk.a],Float64[x0.b,zk.b],color=col)
    if lineas
        PyPlot.plot(tI*[1,1],Float64[x0.a,x0.b],"--",color="black")
        PyPlot.plot(tF*[1,1],Float64[zk.a,zk.b],"--",color="black")
        #PyPlot.xlim(min(0.9*tI,1.1*tI),max(0.9*tF,1.1*tF))
        #PyPlot.ylim(min(0.9*float64(zk.a),1.1*float64(zk.a)),max(0.9*float64(zk.b),1.1*float64(zk.b)))
    end
    return nothing
end
graficaCono(tI,tF,x0::Intervalo,zk::Intervalo)=graficaCono(tI,tF,x0,zk,false,"lightblue")
graficaCono(tI,tF,x0::Intervalo,zk::Intervalo,lineas::Bool)=graficaCono(tI,tF,x0,zk,lineas,"lightblue")
graficaCono(tI,tF,x0::Intervalo,zk::Intervalo,col::ASCIIString)=graficaCono(tI,tF,x0,zk,false,col)
#function graficaCono(tI,tF,x0::Intervalo,zk::Intervalo,col::ASCIIString)
#    PyPlot.fill_between([tI,tF],Float64[x0.a,zk.a],Float64[x0.b,zk.b],color=col)
#    PyPlot.plot(tI*[1,1],Float64[x0.a,x0.b],"--",color="black")
#    PyPlot.plot(tF*[1,1],Float64[zk.a,zk.b],"--",color="black")
#    #PyPlot.xlim(min(0.9*tI,1.1*tI),max(0.9*tF,1.1*tF))
#    #PyPlot.ylim(min(0.9*float64(zk.a),1.1*float64(zk.a)),max(0.9*float64(zk.b),1.1*float64(zk.b)))
#    return nothing
#end
#graficaCono(tI,tF,x0::Intervalo,zk::Intervalo)=graficaCono(tI,tF,x0,zk,"lightblue")







##########################PARTICIONAR EN EL TIEMPO

function envolvente(tI,tF,x0,z0::Intervalo,f::Function,n::Int64)
    # Esta función hace el cálculo para una z0 dada y partiendo el tiempo en n pedazos    
    envolv=Intervalo[]
    tiempo=linspace(tI,tF,2+n)
    x1=zetasKas(tiempo[1],tiempo[2],x0,z0,f)
    push!(envolv,x1)
    for i=1:n
        push!(envolv,conosInterv(tiempo[2+n-1],tiempo[2+n],envolv[i],f))
    end
    return envolv
end

function graficaEnvolvente(tI,tF,x0,envolv::Array{Intervalo,1},col::ASCIIString)
    # Esta función grafica lo que hace la función envolvente
    n=length(envolv)-1
    tiempo=linspace(tI,tF,2+n)
    graficaCono(tiempo[1],tiempo[2],x0,envolv[1],col)
    for i=1:n
        graficaCono(tiempo[1+i],tiempo[2+i],envolv[i],envolv[1+i],col)
    end
    return nothing
end

graficaEnvolvente(tI,tF,x0,envolv::Array{Intervalo,1})=graficaEnvolvente(tI,tF,x0,envolv,"lightblue")




#########EL MÁXIMO EN TODA UNA TRAYECTORIA




function ptosenvolventes(t0,h,n,x0,f,m) 
    T=BuscarTBuena(x0,h,n,f); 
    z0= zetacero(x0,h,n,f); 
    nuevaListaXs=envolvente(t0,T,x0,z0,f,m);
    
    arriba=maximum([(nuevaListaXs[i].b) for i in 1:length(nuevaListaXs)]);
    abajo=minimum([(nuevaListaXs[i].a) for i in 1:length(nuevaListaXs)]);

   return Intervalo(abajo,arriba)
end


function ptosenvolventes(z) 
    arriba=maximum([(z[i].b) for i in 1:length(z)]);
    abajo=minimum([(z[i].a) for i in 1:length(z)]);
   return Intervalo(abajo,arriba)
end

########Resumido para resolver:

function NEuler1Paso(f,x0,t0,n,dt,ordenK)
    tf=t0+dt*n;
    #x_RK=RungeKutta4(t0,x0,dt,n,f);
    T=BuscarTBuena(x0,dt,n,f);
    z0=zetacero(x0,dt,n,f);
    z_k0=kaesimaZeta(t0,T,x0,z0,f,ordenK)[end];
    #z_k=ptosenvolventes(t0,h,n,x0,f,m);
    #graficaEnvolvente(t0,T,x0,nuevaListaXs)
    #PyPlot.plot(linspace(t0,tf,n+1),x_RK,"red");
    #PyPlot.plot(linspace(n*dt,n*dt,n+1),linspace(minimum(x_RK),maximum(x_RK),n+1),"--r")

    nuevaListaXs=envolvente(t0,T,x0,z_k0,f,ordenK)
    return(nuevaListaXs,T)

end





function NEulerMPasos(f,x0,t0,n,dt,ordenK,M)
    cero=[Intervalo(0)];
    xinicial=similar(cero,M);
    tinicial=zeros(M);
    xinicial[1]=Intervalo(x0);
    #tinicial[1]=t0;
    tf=t0+dt*n;

    #primer paso de integración
    nuevaListaXs1=NEuler1Paso(f,x0,t0,n,dt,ordenK)[1];
    T=NEuler1Paso(f,x0,t0,n,dt,ordenK)[2];
    tinicial[1]=T
    z0=zetacero(x0,dt,n,f);

    #segundo paso de integración
    #ya no se usa NEuler1Paso porque se tiene que usar BuscarTBuena partiendo de 
    Xinicial=nuevaListaXs1[end]; #variables mudas
    Tinicial=BuscarTBuena(Xinicial,0,dt,n,f,z0); #variables mudas
    z_k=kaesimaZeta(0,Tinicial,Xinicial,z0,f,ordenK)[end];
    nuevaListaXs2=envolvente(0,Tinicial,Xinicial,z_k,f,ordenK);
    tinicial[2]=tinicial[1]+Tinicial
    xinicial[2]=Xinicial

    # Definiendo el arreglo que aloja las listas de soluciones entre tinicial[i] y 
    #tinicial[i+1]
    nuevaListaXstotal=similar(Array[nuevaListaXs1],M)
    nuevaListaXstotal[1]=nuevaListaXs1
    nuevaListaXstotal[2]=nuevaListaXs2

    #Del tercer paso en adelante
    for i in 3:M
        a=nuevaListaXstotal[i-1]
        Xinicial=a[end];
        Tinicial=BuscarTBuena(Xinicial,0,dt,n,f,z0); #variables mudas
        z_k=kaesimaZeta(0,Tinicial,Xinicial,z0,f,ordenK)[end];
        nuevaListaXstotal[i]=envolvente(0,Tinicial,Xinicial,z_k,f,ordenK);
        xinicial[i]=Xinicial

        tinicial[i]=tinicial[i-1]+Tinicial
        
        #println(Tinicial)
    
        #    tinicial[i]=BuscarTBuena(xinicial[1],0,dt,n,f,z0);
    #z_k=kaesimaZeta(0,tinicial[1],xinicial[1],z0,f,ordenK)[end];
    #nuevaListaXs1=envolvente(0,tinicial[1],xinicial[1],z_k,f,ordenK);
    end

    return(Array[tinicial,xinicial,nuevaListaXstotal])
end


function NEulerMPasos(f,y0::Intervalo,t0,n,dt,ordenK,M)   
    
    cero=[Intervalo(0)];
    xinicial=similar(cero,M);
    tinicial=zeros(M);
    tf=t0+dt*n;
    ###################################bajo############################
    x0=y0.a
    xinicial[1]=Intervalo(x0);
    #tinicial[1]=t0;

    #primer paso de integración
    nuevaListaXs1_bajo=NEuler1Paso(f,x0,t0,n,dt,ordenK)[1];
    T=NEuler1Paso(f,x0,t0,n,dt,ordenK)[2];
    tinicial_bajo=T
    z0_bajo=zetacero(x0,dt,n,f);
    
    ##################################alto##############################
    x0=y0.b
    xinicial[1]=Intervalo(x0);
    #tinicial[1]=t0;

    #primer paso de integración
    nuevaListaXs1_alto=NEuler1Paso(f,x0,t0,n,dt,ordenK)[1];
    T=NEuler1Paso(f,x0,t0,n,dt,ordenK)[2];
    tinicial_alto=T
    z0_alto=zetacero(x0,dt,n,f);

    z0=uneInterv([z0_bajo],[z0_alto])[1];
    tinicial[1]=minimum([tinicial_bajo,tinicial_alto])
    if tinicial[1]==tinicial_bajo
        nuevaListaXs1=nuevaListaXs1_bajo
    else
        nuevaListaXs1=nuevaListaXs1_alto
    end



    #segundo paso de integración
    #ya no se usa NEuler1Paso porque se tiene que usar BuscarTBuena partiendo de 
    Xinicial=nuevaListaXs1[end]; #variables mudas
    Tinicial=BuscarTBuena(Xinicial,0,dt,n,f,z0); #variables mudas
    z_k=kaesimaZeta(0,Tinicial,Xinicial,z0,f,ordenK)[end];
    nuevaListaXs2=envolvente(0,Tinicial,Xinicial,z_k,f,ordenK);
    tinicial[2]=tinicial[1]+Tinicial
    xinicial[2]=Xinicial

    # Definiendo el arreglo que aloja las listas de soluciones entre tinicial[i] y 
    #tinicial[i+1]
    nuevaListaXstotal=similar(Array[nuevaListaXs1],M)
    nuevaListaXstotal[1]=nuevaListaXs1
    nuevaListaXstotal[2]=nuevaListaXs2

    #Del tercer paso en adelante
    for i in 3:M
        a=nuevaListaXstotal[i-1]
        Xinicial=a[end];
        Tinicial=BuscarTBuena(Xinicial,0,dt,n,f,z0); #variables mudas
        z_k=kaesimaZeta(0,Tinicial,Xinicial,z0,f,ordenK)[end];
        nuevaListaXstotal[i]=envolvente(0,Tinicial,Xinicial,z_k,f,ordenK);
        xinicial[i]=Xinicial

        tinicial[i]=tinicial[i-1]+Tinicial
        
        #println(Tinicial)
    
        #    tinicial[i]=BuscarTBuena(xinicial[1],0,dt,n,f,z0);
    #z_k=kaesimaZeta(0,tinicial[1],xinicial[1],z0,f,ordenK)[end];
    #nuevaListaXs1=envolvente(0,tinicial[1],xinicial[1],z_k,f,ordenK);
    end

    return(Array[tinicial,xinicial,nuevaListaXstotal])
end




########## Colores para graficar

cnames = [
"aliceblue",            "#F0F8FF",
"antiquewhite",         "#FAEBD7",
"aqua",                 "#00FFFF",
"aquamarine",           "#7FFFD4",
"azure",                "#F0FFFF",
"beige",                "#F5F5DC",
"bisque",               "#FFE4C4",
"black",                "#000000",
"blanchedalmond",       "#FFEBCD",
"blue",                 "#0000FF",
"blueviolet",           "#8A2BE2",
"brown",                "#A52A2A",
"burlywood",            "#DEB887",
"cadetblue",            "#5F9EA0",
"chartreuse",           "#7FFF00",
"chocolate",            "#D2691E",
"coral",                "#FF7F50",
"cornflowerblue",       "#6495ED",
"cornsilk",             "#FFF8DC",
"crimson",              "#DC143C",
"cyan",                 "#00FFFF",
"darkblue",             "#00008B",
"darkcyan",             "#008B8B",
"darkgoldenrod",        "#B8860B",
"darkgray",             "#A9A9A9",
"darkgreen",            "#006400",
"darkkhaki",            "#BDB76B",
"darkmagenta",          "#8B008B",
"darkolivegreen",       "#556B2F",
"darkorange",           "#FF8C00",
"darkorchid",           "#9932CC",
"darkred",              "#8B0000",
"darksalmon",           "#E9967A",
"darkseagreen",         "#8FBC8F",
"darkslateblue",        "#483D8B",
"darkslategray",        "#2F4F4F",
"darkturquoise",        "#00CED1",
"darkviolet",           "#9400D3",
"deeppink",             "#FF1493",
"deepskyblue",          "#00BFFF",
"dimgray",              "#696969",
"dodgerblue",           "#1E90FF",
"firebrick",            "#B22222",
"floralwhite",          "#FFFAF0",
"forestgreen",          "#228B22",
"fuchsia",              "#FF00FF",
"gainsboro",            "#DCDCDC",
"ghostwhite",           "#F8F8FF",
"gold",                 "#FFD700",
"goldenrod",            "#DAA520",
"gray",                 "#808080",
"green",                "#008000",
"greenyellow",          "#ADFF2F",
"honeydew",             "#F0FFF0",
"hotpink",              "#FF69B4",
"indianred",            "#CD5C5C",
"indigo",               "#4B0082",
"ivory",                "#FFFFF0",
"khaki",                "#F0E68C",
"lavender",             "#E6E6FA",
"lavenderblush",        "#FFF0F5",
"lawngreen",            "#7CFC00",
"lemonchiffon",         "#FFFACD",
"lightblue",            "#ADD8E6",
"lightcoral",           "#F08080",
"lightcyan",            "#E0FFFF",
"lightgoldenrodyellow", "#FAFAD2",
"lightgreen",           "#90EE90",
"lightgray",            "#D3D3D3",
"lightpink",            "#FFB6C1",
"lightsalmon",          "#FFA07A",
"lightseagreen",        "#20B2AA",
"lightskyblue",         "#87CEFA",
"lightslategray",       "#778899",
"lightsteelblue",       "#B0C4DE",
"lightyellow",          "#FFFFE0",
"lime",                 "#00FF00",
"limegreen",            "#32CD32",
"linen",                "#FAF0E6",
"magenta",              "#FF00FF",
"maroon",               "#800000",
"mediumaquamarine",     "#66CDAA",
"mediumblue",           "#0000CD",
"mediumorchid",         "#BA55D3",
"mediumpurple",         "#9370DB",
"mediumseagreen",       "#3CB371",
"mediumslateblue",      "#7B68EE",
"mediumspringgreen",    "#00FA9A",
"mediumturquoise",      "#48D1CC",
"mediumvioletred",      "#C71585",
"midnightblue",         "#191970",
"mintcream",            "#F5FFFA",
"mistyrose",            "#FFE4E1",
"moccasin",             "#FFE4B5",
"navajowhite",          "#FFDEAD",
"navy",                 "#000080",
"oldlace",              "#FDF5E6",
"olive",                "#808000",
"olivedrab",            "#6B8E23",
"orange",               "#FFA500",
"orangered",            "#FF4500",
"orchid",               "#DA70D6",
"palegoldenrod",        "#EEE8AA",
"palegreen",            "#98FB98",
"paleturquoise",        "#AFEEEE",
"palevioletred",        "#DB7093",
"papayawhip",           "#FFEFD5",
"peachpuff",            "#FFDAB9",
"peru",                 "#CD853F",
"pink",                 "#FFC0CB",
"plum",                 "#DDA0DD",
"powderblue",           "#B0E0E6",
"purple",               "#800080",
"red",                  "#FF0000",
"rosybrown",            "#BC8F8F",
"royalblue",            "#4169E1",
"saddlebrown",          "#8B4513",
"salmon",               "#FA8072",
"sandybrown",           "#FAA460",
"seagreen",             "#2E8B57",
"seashell",             "#FFF5EE",
"sienna",               "#A0522D",
"silver",               "#C0C0C0",
"skyblue",              "#87CEEB",
"slateblue",            "#6A5ACD",
"slategray",            "#708090",
"snow",                 "#FFFAFA",
"springgreen",          "#00FF7F",
"steelblue",            "#4682B4",
"tan",                  "#D2B48C",
"teal",                 "#008080",
"thistle",              "#D8BFD8",
"tomato",               "#FF6347",
"turquoise",            "#40E0D0",
"violet",               "#EE82EE",
"wheat",                "#F5DEB3",
"white",                "#FFFFFF",
"whitesmoke",           "#F5F5F5",
"yellow",               "#FFFF00",
"yellowgreen",          "#9ACD32"];

colores=[cnames[2*i+1] for i in 0:length(cnames)/2-1];





######UNIDO PARA GRAFICAR

#function graficarEulerM(tinicial,xinicial,nuevaListaXstotal,x_RK,t0,tf,n,dt,pasosdeintegracion)
#    PyPlot.plot(linspace(t0,tf,n+1),x_RK,"red");
#    PyPlot.plot(linspace(n*dt,n*dt,n+1),linspace(minimum(x_RK),maximum(x_RK),n+1),"--r")
#
#    graficaEnvolvente(t0,tinicial[1],xinicial[1],nuevaListaXstotal[1],colores[1])
#    graficaEnvolvente(tinicial[1],tinicial[2],xinicial[2],nuevaListaXstotal[2],colores[2])
#
#    for i in 3:pasosdeintegracion
#        graficaEnvolvente(tinicial[i-1],tinicial[i],xinicial[i],nuevaListaXstotal[i], colores[i])
#    end
#    println("el número máximo alcanzado con intervalos es: $(tinicial[end]) de $(tf)")    
#end

#function graficarEulerM(tinicial,xinicial,nuevaListaXstotal,x_RK,t0,tf,n,dt,pasosdeintegracion,graficartodoRungeKutta::Bool)
#
#    if graficartodoRungeKutta==true
#        graficarEulerM(tinicial,xinicial,nuevaListaXstotal,x_RK,t0,tf,n,dt,pasosdeintegracion)
#
#    else 
#        PyPlot.plot(linspace(t0,tf,n+1),x_RK,"red");
#
#        graficaEnvolvente(t0,tinicial[1],xinicial[1],nuevaListaXstotal[1],colores[1])
#        graficaEnvolvente(tinicial[1],tinicial[2],xinicial[2],nuevaListaXstotal[2],colores[2])
#       
#        println("el número máximo alcanzado con intervalos es: $(tinicial[end]) de $(tf)")    
#    
#       
#            for i in 3:pasosdeintegracion
#                graficaEnvolvente(tinicial[i-1],tinicial[i],xinicial[i],nuevaListaXstotal[i], colores[i])    
#            end
#
#                println("Proporcione un número de pasos (de longitud $(dt)) extra para graficar Runge Kuta después del método con intervalos")    
#                n_extra=readline()
#                n_extra=int(n_extra)
#                n_extra + length(tinicial) > n && error("número de pasos extras demasiado grande")
#
#                PyPlot.xlim([t0,tinicial[end]+n_extra*dt]); 
#    end     
#end
    

function graficarEulerM(tinicial,xinicial,nuevaListaXstotal,x_RK,t0,tf,n,dt,pasosdeintegracion,graficartodoRungeKutta::Bool)
    pasosdeintegracion>140 && error("Te pasaste de colores")
    PyPlot.plot(linspace(t0,tf,n+1),x_RK,"red");
    graficaEnvolvente(t0,tinicial[1],xinicial[1],nuevaListaXstotal[1],colores[1])
    graficaEnvolvente(tinicial[1],tinicial[2],xinicial[2],nuevaListaXstotal[2],colores[2])
    println("el número máximo alcanzado con intervalos es: $(tinicial[end]) de $(tf)")    
    for i in 3:pasosdeintegracion
        graficaEnvolvente(tinicial[i-1],tinicial[i],xinicial[i],nuevaListaXstotal[i], colores[i])    
    end
    if graficartodoRungeKutta==false
        println("Proporcione un número de pasos (de longitud $(dt)) extra para graficar Runge Kuta después del método con intervalos")    
        n_extra=readline()
        n_extra=int(n_extra)
        n_extra + length(tinicial) > n && error("número de pasos extras demasiado grande")
        PyPlot.xlim([t0,tinicial[end]+n_extra*dt]); 
    end     
end

graficarEulerM(tinicial,xinicial,nuevaListaXstotal,x_RK,t0,tf,n,dt,pasosdeintegracion)=graficarEulerM(tinicial,xinicial,nuevaListaXstotal,x_RK,t0,tf,n,dt,pasosdeintegracion,true)



end

