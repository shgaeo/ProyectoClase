###

module AutoDif

using Intervalos

import Base.+,Base.-,Base.*,Base./,Base.^,Base.==,Base.in
import Base.exp, Base.log, Base.sin, Base.cos, Base.tan

export ValorDeriv


type ValorDeriv
    valor
    deriv
end

ValorDeriv(valor)= ValorDeriv(valor,1.0)
ValorDeriv(x::Intervalo)=ValorDeriv(x,Intervalo(1.0))

+(x::ValorDeriv,y::ValorDeriv)=ValorDeriv(x.valor+y.valor , x.deriv+y.deriv) #defino la suma entre ValorDeriv
+(x::Number,y::ValorDeriv)=ValorDeriv(x,0.0)+y #defino la suma número + ValorDeriv
+(y::ValorDeriv,x::Number)=x+y #defino la suma ValorDeriv + número

-(x::ValorDeriv,y::ValorDeriv)=ValorDeriv(x.valor-y.valor , x.deriv-y.deriv) #defino la resta entre ValorDeriv
-(x::Number,y::ValorDeriv)=ValorDeriv(x,0.0)-y #defino la resta número - ValorDeriv
-(y::ValorDeriv,x::Number)=y-ValorDeriv(x,-0.0) #defino la suma ValorDeriv - número

*(x::ValorDeriv,y::ValorDeriv)=ValorDeriv(x.valor*y.valor , x.valor*y.deriv+x.deriv*y.valor) #defino producto entre ValorDeriv
*(x::Number,y::ValorDeriv)=ValorDeriv(x,0.0)*y #defino producto ValorDeriv por número
*(y::ValorDeriv,x::Number)=x*y #defino producto ValorDeriv por número

/(x::ValorDeriv,y::ValorDeriv)=ValorDeriv(x.valor/y.valor , (x.deriv*y.valor-x.valor*y.deriv)/y.valor^2) #defino división
/(x::Number,y::ValorDeriv)=ValorDeriv(x,0.0)/y
/(x::ValorDeriv,y::Number)=x/ValorDeriv(y,0.0)

^(x::ValorDeriv,y::Integer)=ValorDeriv(x.valor^y , y*x.valor^(y-1)*x.deriv) #defino potencia
^(x::ValorDeriv,y::Number)=ValorDeriv(x.valor^y , y*x.valor^(y-1)*x.deriv) #defino potencia

exp(x::ValorDeriv)=ValorDeriv(exp(x.valor) , exp(x.valor)*x.deriv) #la derivada de exp es ella misma

log(x::ValorDeriv)=ValorDeriv(log(x.valor) , (1/(x.valor))*x.deriv) #la derivada de log(x) es 1/x

sin(x::ValorDeriv)=ValorDeriv(sin(x.valor) , cos(x.valor)*x.deriv) #la derivada de sin(x) es cos(x)

cos(x::ValorDeriv)=ValorDeriv(cos(x.valor) , -sin(x.valor)*x.deriv) #la derivada de cos(x) es -sin(x)

tan(x::ValorDeriv)=ValorDeriv(tan(x.valor) , (1/(cos(x.valor))^2)*x.deriv) #la derivada de tan(x) es 1/(cos(x))^2

end
