###

module Intervalos

import Base.show
import Base.+,Base.-,Base.*,Base./,Base.^,Base.==,Base.in
import Base.exp, Base.log, Base.atan, Base.asin, Base.acos

export Intervalo, redonUP, redonDOWN, uNir, iNtersectar

type Intervalo
    a::BigFloat
    b::BigFloat
    
    Intervalo(a::BigFloat,b::BigFloat) = ( a>=b ? new(b,a) : new(a,b))
end

function redonUP(x)
    with_rounding(BigFloat,RoundUp) do 
        BigFloat(x) 
    end
end

function redonDOWN(x)
    with_rounding(BigFloat,RoundDown) do 
        BigFloat(x) 
    end
end

Intervalo(a,b)=Intervalo(redonDOWN(a),redonUP(b)) #Con esta definición puedo usar como argumento cadenas de caracteres
                                                    # tomando en cuenta el redondeo
                                                    

Intervalo(a)=Intervalo(a,a)

show(io::IO, a::Intervalo) = print(io::IO, "[$(a.a),
$(a.b)]")

function redonUP(f::Function,x,y)
    with_rounding(BigFloat,RoundUp) do 
        f(BigFloat(x),BigFloat(y))
    end
end

function redonDOWN(f::Function,x,y)
    with_rounding(BigFloat,RoundDown) do 
        f(BigFloat(x),BigFloat(y))
    end
end

function redonUP(f::Function,x)
    with_rounding(BigFloat,RoundUp) do 
        f(BigFloat(x))
    end
end

function redonDOWN(f::Function,x)
    with_rounding(BigFloat,RoundDown) do 
        f(BigFloat(x))
    end
end

function definidor!(f::Function)  #función que define operaciones entre intervalos
    if f!=^
        f(x::Intervalo,y::Intervalo) = (f==/ && y.a<0.0 && y.b>0.0) ? error("El intervalo divisor no puede contener al cero") :
            Intervalo(min(redonDOWN(f,x.a,y.a),redonDOWN(f,x.b,y.b),redonDOWN(f,x.a,y.b),redonDOWN(f,x.b,y.a),
                redonUP(f,x.a,y.a),redonUP(f,x.b,y.b),redonUP(f,x.a,y.b),redonUP(f,x.b,y.a)),
            max(redonDOWN(f,x.a,y.a),redonDOWN(f,x.b,y.b),redonDOWN(f,x.a,y.b),redonDOWN(f,x.b,y.a),
                redonUP(f,x.a,y.a),redonUP(f,x.b,y.b),redonUP(f,x.a,y.b),redonUP(f,x.b,y.a)))
        f(x::Intervalo,y::Number) = (f==/ && (redonDOWN(y)==0.0 || redonUP(y)==0.0)) ?
            error("No puedo dividir un intervalo entre cero") :
            f(x,Intervalo(y))
        f(x::Intervalo,y::ASCIIString) = (f==/ && (redonDOWN(y)==0.0 || redonUP(y)==0.0)) ?
            error("No puedo dividir un intervalo entre cero") :
            f(x,Intervalo(y))
        f(y::Number,x::Intervalo) = (f==/ && x.a<0.0 && x.b>0.0) ? error("El intervalo divisor no puede contener al cero") :
            f(Intervalo(y),x)
        f(y::ASCIIString,x::Intervalo) = (f==/ && x.a<0.0 && x.b>0.0) ? error("El intervalo divisor no puede contener al cero") :
            f(Intervalo(y),x)
     else
        f(x::Intervalo,y::Integer)= y<0 ? "1"/f(x,-y) :
            (y==0 ? Intervalo("1") :
            (x.a<0.0 && x.b>0.0 ? Intervalo(min(0.0,redonDOWN(f,x.a,y),redonDOWN(f,x.b,y)),max(redonDOWN(f,x.a,y),redonDOWN(f,x.b,y))) :
            Intervalo(min(redonDOWN(f,x.a,y),redonDOWN(f,x.b,y)),max(redonDOWN(f,x.a,y),redonDOWN(f,x.b,y)))))       
    end
end

definidor!(+)
definidor!(-)
-(x::Intervalo) = x * -1
definidor!(*)
definidor!(/)
definidor!(^)

==(x::Intervalo,y::Intervalo) = (x.a==y.a && x.b==y.b) ? true : false
in(x,y::Intervalo) = (y.a<=Intervalo(x).a && Intervalo(x).b<=y.b) ? true : false
in(x::Intervalo,y::Intervalo) = in(x.a,y) && in(x.b,y)
exp(x::Intervalo)=Intervalo(redonDOWN(exp,x.a),redonUP(exp,x.b))
log(x::Intervalo)= x.a<0 ? error("Función log no acepta números negativos \n") : Intervalo(redonDOWN(log,x.a),redonUP(log,x.b))
atan(x::Intervalo)=Intervalo(redonDOWN(atan,x.a),redonUP(atan,x.b))
asin(x::Intervalo)=
        x.a<-1||x.b>1 ? error("Función asin solo acepta números entre -1 y 1 \n") : Intervalo(redonDOWN(asin,x.a),redonUP(asin,x.b))
acos(x::Intervalo)=
        x.a<-1||x.b>1 ? error("Función acos solo acepta números entre -1 y 1 \n") : Intervalo(redonDOWN(acos,x.b),redonUP(acos,x.a))
^(x::Intervalo,y)=exp(y*log(x))

uNir(x::Intervalo,y::Intervalo) = (x.a<=y.b && y.a<=x.b) ? Intervalo(min(x.a,y.a),max(x.b,y.b)) : return("Deben ser intervalos que se intersectan")
iNtersectar(x::Intervalo,y::Intervalo) = (x.a<=y.b && y.a<=x.b) ? Intervalo(max(x.a,y.a),min(x.b,y.b)) : return("Deben ser intervalos que se intersectan")


end
