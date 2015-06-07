###

module Graficador

export graficador, graficaarea

using PyPlot, Intervalos    

function graficador(a,b,f::Function)
    puntosX=(-5:0.125:8)
    puntosY=similar(puntosX)
    for i=1:length(puntosX)
        puntosY[i]=f(puntosX[i])
    end
    
    ptsX=(a:0.125:b)
    ptsY=similar(ptsX)
    for i=1:length(ptsX)
        ptsY[i]=f(ptsX[i])
    end
    
    polY=f(Intervalo(a,b))
    
    PyPlot.plot(puntosX,puntosY,"r",linewidth=1)
    PyPlot.plot(ptsX, ptsY, "b", linewidth=2.0)
    
    PyPlot.plot((a-1)*[1,1],float64([polY.a,polY.b]), "black", linewidth=2.0)
    PyPlot.plot([a-1,a],float64(polY.a)*[1,1],"black", linewidth=1.0, linestyle="--")
    PyPlot.plot([a-1,a],float64(polY.b)*[1,1], "black", linewidth=1.0, linestyle="--")
    
    PyPlot.plot([a,b],(float64(polY.a)-1)*[1,1], "black", linewidth=2.0)
    PyPlot.plot([a,a],[float64(polY.a)-1,float64(polY.a)], "black", linewidth=1.0, linestyle="--")
    PyPlot.plot([b,b],[float64(polY.a)-1,float64(polY.a)], "black", linewidth=1.0, linestyle="--")
    
    PyPlot.xlim(min(-1,a-1),max(0,b+1))
    PyPlot.ylim(min(-1,float64(polY.a)-1.5),max(float64(polY.a),float64(polY.b),f(a),f(b))+1);

    PyPlot.fill_between([a,b],float64(polY.a),float64(polY.b),color="lightgray")
    PyPlot.suptitle("La sombra gris es la zona calculada con intervalo [$a,$b]")
end

function graficaarea(a,b,f::Function)
    puntosX=(a-2:0.125:b+2)
    puntosY=similar(puntosX)
    for i=1:length(puntosX)
        puntosY[i]=f(puntosX[i])
    end
    
    ptsX=(a:0.125:b)
    ptsY=similar(ptsX)
    for i=1:length(ptsX)
        ptsY[i]=f(ptsX[i])
    end
    
    polY=f(Intervalo(a,b))
    
    PyPlot.plot(puntosX,puntosY,"r",linewidth=1)
    PyPlot.plot(ptsX, ptsY, "b", linewidth=2.0)
        
    PyPlot.xlim(min(-1,a-1),max(0,b+1))
    PyPlot.ylim(min(-1,float64(polY.a)-1.5),max(float64(polY.a),float64(polY.b),f(a),f(b))+1);

    PyPlot.fill_between([a,b],float64(polY.a),float64(polY.b),color="lightgray")
end

end
