from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from pymatgen.core import Composition
from typing import Union, List
from importlib import resources

base = importr('base')
utils = importr('utils')
locfit = importr('locfit')

r = robjects.r
path = str(resources.files('pqam_dparamhu2021'))
r['source'](path+'/HEA_pred.R')
# Initialize the models
heaPredInit = robjects.globalenv['init']
heaPredInit(path)

# Load the prediction function
heaPredFunc = robjects.globalenv['HEA_pred']

# (Ti,Zr,Hf,V,Nb,Ta,Mo,W,Re,Ru)
elementsSpace = ['Ti', 'Zr', 'Hf', 'V', 'Nb', 'Ta', 'Mo', 'W', 'Re', 'Ru']


def predict(comp: Union[str, Composition]) -> list:
    """
    Predicts the GSF, Surd, and resulting D parameter for a given HEA composition in the
    composition space of (Ti,Zr,Hf,V,Nb,Ta,Mo,W,Re,Ru) based on Hu's 2021 model (10.1016/j.actamat.2021.116800).

    Args:
        comp: A composition string which will be cast into pymatgen Composition object or ready Composition object.

    Returns:
        A float list representing the predicted GSF, Surd, and D parameter.
    """

    assert isinstance(comp, (str, Composition)), \
        "comp must be a string or a pymatgen Composition object."
    if isinstance(comp, str):
        comp = Composition(comp)

    assert all([e.symbol in elementsSpace for e in comp.elements]), \
        "The composition must be in the composition space of (Ti,Zr,Hf,V,Nb,Ta,Mo,W,Re,Ru)."

    compList = [comp.get_atomic_fraction(e) for e in elementsSpace]
    result = heaPredFunc(robjects.FloatVector(compList), path)
    return list(result)


def cite() -> List[str]:
    """
    Returns the citation of the model.

    Returns:
        A list of strings representing the citation of the model.
    """
    return ["Yong-Jie Hu, Aditya Sundar, Shigenobu Ogata, Liang Qi, Screening of generalized stacking fault energies, "
            "surface energies and intrinsic ductile potency of refractory multicomponent alloys, Acta Materialia, "
            "Volume 210, 2021, 116800, https://doi.org/10.1016/j.actamat.2021.116800."
            ]
