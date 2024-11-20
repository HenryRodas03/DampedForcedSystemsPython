from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import sympy as sp

app = FastAPI()

# Configurar CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:4200"],  # Agrega el origen de tu app Angular
    allow_credentials=True,
    allow_methods=["*"],  # Permitir todos los métodos (GET, POST, etc.)
    allow_headers=["*"],  # Permitir todos los encabezados
)

# Modelo para la entrada de datos
class SimulationInput(BaseModel):
    m: float
    k: float
    b: float
    c1: float
    c2: float
    ftType: int
    fT1: float
    fT2: float

def particularFunction (m: float, k: float, b: float,ftType: int, fT1: float, fT2: float):
    t = sp.symbols('t')  # Variable de tiempo

    # Calcular valores
    if ftType == 1:

        ft = 0
        return ft

    elif ftType ==2:
        # Definir fuerza externa f(t) usando fT1 y fT2
        ft = fT1 * sp.sin(fT2 * t)

        # Resolver la solución particular de la ecuación diferencial
        A, B = sp.symbols('A B', real=True)  # Coeficientes a determinar
        yp = A * sp.sin(fT2 * t) + B * sp.cos(fT2 * t)  # Forma tentativa de la solución particular

        # Derivadas
        yp_prime = sp.diff(yp, t)
        yp_double_prime = sp.diff(yp_prime, t)

        # Sustitución en la ecuación diferencial
        lhs = yp_double_prime + (b/m) * yp_prime + (k/m) * yp - (ft/m)

        # Resolver para A y B
        equations = [lhs.coeff(sp.sin(fT2 * t)), lhs.coeff(sp.cos(fT2 * t))]
        coefficients = sp.solve(equations, (A, B))

        # Redondear los coeficientes a 3 decimales
        A_value = round(float(coefficients[A]), 4)
        B_value = round(float(coefficients[B]), 4)

        # Construir la solución particular y simplificar
        particular_solution = A_value * sp.sin(fT2 * t) + B_value * sp.cos(fT2 * t)

        # Retornar la solución total como cadena
        return particular_solution

    else:

        ft = t

        A, B = sp.symbols('A B', real=True)  # Coeficientes a determinar
        yp = A * t + B   # Forma tentativa de la solución particular
        
        # Derivadas
        yp_prime = sp.diff(yp, t)
        yp_double_prime = sp.diff(yp_prime, t)

        # Sustitución en la ecuación diferencial
        lhs = yp_double_prime + (b/m) * yp_prime + (k/m) * yp - (ft/m)

        # Resolver para A y B
        equations = sp.collect(lhs, [A, B])  # Agrupar términos dependientes de A y B
        coefficients = sp.solve(equations, (A, B))

        # Verificar si se encontró una solución
        if isinstance(coefficients, dict):
            # Si sp.solve devuelve un diccionario
            A_value = round(float(coefficients[A]), 4)
            B_value = round(float(coefficients[B]), 4)
        elif isinstance(coefficients, list) and len(coefficients) > 0:
            # Si sp.solve devuelve una lista, tomamos la primera solución
            A_value = round(float(coefficients[0][0]), 4)
            B_value = round(float(coefficients[0][1]), 4)
        else:
            raise ValueError("No se encontró solución para los coeficientes A y B")

        # Construir la solución particular
        particular_solution = A_value * t + B_value

        # Retornar la solución particular
        return particular_solution



def solve_motion(m: float, k: float, b: float, c1: float, c2: float, ftType: int, fT1: float, fT2: float):
    t = sp.symbols('t')  # Variable de tiempo
    sistemtype = ""
    # Convertir fT1 y fT2 a float si no lo son
    fT1 = float(fT1)
    fT2 = float(fT2)

    # Cálculo de la expresión para resolv
    resolv = ((b / (2 * m))**2) - (k / m)

    # Condición para modificar el valor de resolv
    if resolv > 0:

        l1 = round(float((-b / (2 * m)) + sp.sqrt(resolv)), 3)
        l2 = round(float((-b / (2 * m)) - sp.sqrt(resolv)), 3)

        # Construir la solución complementaria
        complementary_solution = c1 * sp.exp(l1 * t) + c2 * sp.exp(l2 * t)

        vt = sp.diff(complementary_solution, t)
        at = sp.diff(vt, t)

        particular_solution = particularFunction(m,k,b,ftType,fT1,fT2)

        total_solution = complementary_solution + particular_solution
        
        sistemtype = "Sobreamortiguado"

    elif resolv == 0:
        l = -b / (2 * m)

        # Construir la solución complementaria
        complementary_solution = c1 * t * sp.exp(l * t) + c2 * sp.exp(l * t)

        particular_solution = particularFunction(m,k,b,ftType,fT1,fT2)

        vt = sp.diff(complementary_solution, t)
        at = sp.diff(vt, t)

        total_solution = complementary_solution + particular_solution

        sistemtype = "Criticamenteamortiguado"

        
    else:
        a = round(b / (2 * m), 3)

        wd = round(float(sp.sqrt(abs(((k/m)**2)-(a**2)))), 3)
  
        # Construir la solución complementaria
        complementary_solution = sp.exp(-a*t) * (c1 * sp.cos(wd * t) + c2 * sp.sin(wd * t))

        vt = sp.diff(complementary_solution, t)
        at = sp.diff(vt, t)

        particular_solution = particularFunction(m,k,b,ftType,fT1,fT2)

        total_solution = complementary_solution + particular_solution

        sistemtype = "Subamortiguado"

    # Retornar la solución total como cadena
    return {
                "status": "OK",
                "data": {
                "sistemtype": sistemtype,
                    "yt": str(total_solution),  # Convertir a cadena para hacerla serializable
                    "vt":str(vt),
                    "at":str(at)
                }
            }        

@app.post("/simulate")
async def simulate(input_data: SimulationInput):
    result = solve_motion(
        input_data.m,
        input_data.k,
        input_data.b,
        input_data.c1,
        input_data.c2,
        input_data.ftType,
        input_data.fT1,
        input_data.fT2,
    )
    return result
