
    module types
    
    
    
    type planet
        real:: r(2), v(2), m
    end type
        
    type system
        type(planet):: bodies(3)
    end type
    
    type gravity_planet
        real:: drdt(2), dvdt(2)
        real:: ai(2)
    end type
    
    end module