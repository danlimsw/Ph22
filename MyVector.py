class MyVector(list):
    """A list-based vector class"""

    def __add__(self,other):
        """Element-by-element addition, or addition of constant"""
        try:
            return MyVector(map(lambda x,y:x+y,self,other))
        except TypeError:
            return MyVector(map(lambda x: x+other,self))
    def __neg__(self):
        """Unary negation"""
        try:
            return MyVector(map(lambda x:-x,self))
        except TypeError:
            return MyVector(-self)
    def __sub__(self,other):
        """Element-by-element subtraction, or subtraction of constant"""
        try:
            return MyVector(map(lambda x,y:x-y,self,other))
        except TypeError:
            return MyVector(map(lambda x: x-other,self))
    def __mul__(self,other):
        """Element-by-element multiplication (dot product), or
            multiplication by constant"""
        try:
            return MyVector(map(lambda x,y:x*y,self,other))
        except TypeError:
            return MyVector(map(lambda x:x*other,self))
    def __div__(self,other):
        """Element-by-element division, or division by constant"""
        try:
            return MyVector(map(lambda x,y:float(x)/float(y),self,other))
        except TypeError:
            return MyVector(map(lambda x:float(x)/float(other),self))
    def __radd__(self,other):
        """Right addition by constant"""
        return MyVector(map(lambda x:x+other,self))
    def __rsub__(self,other):
        """Right subtractions by constant"""
        return MyVector(map(lambda x:other-x,self))
    def __rmul__(self,other):
        """Right multiplication by constant"""
        return MyVector(map(lambda x: x*other,self))
    def __rdiv__(self,other):
        """Right division by constant"""
        return MyVector(map(lambda x:float(other)/float(x),self))
    def cross(u,v):
        try:
            return MyVector([u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]])
        except TypeError:
            return 'Wrong dimensions/type'
