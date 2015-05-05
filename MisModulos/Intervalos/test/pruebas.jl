using Intervalos,FactCheck

facts("Probando tipo Intervalo") do
    interv1=Intervalo("0.1")
    @fact typeof(interv1.a) => BigFloat 
    @fact Intervalo("4","3") => Intervalo(3,4)
    @fact -Intervalo("4","3") => Intervalo(-4.0,-3.0)
    @fact interv1 => Intervalo(1//10) #verifico que convierta cadena de caracteres y racionales
    interv2=Intervalo("0.1",10//3)
    numero1=with_rounding(BigFloat,RoundDown) do
            BigFloat("0.1")
        end
    numero2=with_rounding(BigFloat,RoundUp) do
            BigFloat(10//3)
        end
    @fact interv2.a => numero1 #verifico que haga bien el redondeo
    @fact interv2.b => numero2 #verifico que haga bien el redondeo
    @fact Intervalo("0.2").a => less_than(Intervalo("0.2").b) #verifico que haga bien el redondeo
    @fact in("0.1",Intervalo("0.1")) => true
end

facts("Probando operaciones entre intervalos")do
    interv1=Intervalo("-1.0","3.0")
    interv2=Intervalo("12","13" )
    @fact interv1+interv2 => Intervalo("11","16") #verifico operaciones entre intervalos "enteros"
    @fact interv1-interv2 => Intervalo("-14","-9")
    @fact interv1*interv2 => Intervalo("-13","39")
    @fact interv1/interv2 => Intervalo(-1//12,3//12)
    numero1=4
    @fact interv1+numero1 => Intervalo("3","7") #verifico operaciones entre intervalos y números
    @fact interv1-numero1 => Intervalo("-5","-1")
    @fact interv1*numero1 => Intervalo("-4","12")
    @fact interv2/numero1 => Intervalo("3",13//4)
    @fact numero1/interv2 => Intervalo(4//13,4//12)
    @fact interv1^numero1 => Intervalo("0","81")
    @fact interv2^-1 => 1/interv2
    interv3=Intervalo("1","4")
    interv4=Intervalo("0.2","3")
    numero2=with_rounding(BigFloat,RoundDown) do
            BigFloat("1.2")
        end
    numero3=10//3
    numero4=with_rounding(BigFloat,RoundUp) do
            BigFloat(2//3)
        end
    @fact (interv3+interv4).a => less_than_or_equal(numero2) #pruebo redondeo en operaciones
    @fact (interv3-numero3).b => greater_than_or_equal(numero4)    
end

facts("Probando funciones sobre intervalos")do
    interv1=Intervalo("0.0","0.5")
    @fact exp(BigFloat(1.0)) => exp(Intervalo("1.0")).b
    @fact exp(Intervalo(1)).a => less_than(exp(Intervalo(1)).b) #checo redondeo
    @fact exp(interv1).a => 1.0
    @fact exp(interv1).b => greater_than_or_equal(exp(interv1.b))
    
    @fact log(interv1).a => -inf(BigFloat)
    @fact log(interv1).b => greater_than_or_equal(log(interv1.b))
    @fact log(interv1).b => less_than(1.0)
    
    @fact atan(interv1).a => 0.0
    @fact atan(interv1).b => greater_than_or_equal(atan(interv1.b))
    @fact atan(Intervalo(inf(BigFloat))).a => π/BigFloat(2)
    @fact atan(Intervalo(-inf(BigFloat))).b => π/BigFloat(-2)
    
    @fact asin(interv1).a => 0.0
    @fact asin(interv1).b => greater_than_or_equal(asin(interv1.b))
    @fact asin(Intervalo(1)).a => π/BigFloat(2)
    @fact asin(Intervalo(-1)).b => π/BigFloat(-2)
    
    @fact acos(interv1).a => less_than_or_equal(acos(interv1.b))
    @fact acos(interv1).b => nextfloat(π/BigFloat(2))
    @fact acos(Intervalo(1)) => Intervalo(0,0)
    @fact acos(Intervalo(-1)).a => π
    
    @fact (interv1^"0.5").a => 0.0 
    @fact (interv1^3.1).b => greater_than_or_equal((interv1.b)^3.1) #si este manda error, falta incluir redondeo
end;


