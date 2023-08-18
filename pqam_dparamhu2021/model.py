from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from rpy2.robjects import conversion, default_converter
from pymatgen.core import Composition
from typing import Union, List
from importlib import resources
import sys

base = importr('base')
utils = importr('utils')
locfit = importr('locfit')

with conversion.localconverter(default_converter):
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


def predict(
        comp: Union[Composition, str],
        outputType: str = "array") -> Union[dict, list]:
    """
    Predicts the GSF, Surf, and resulting D parameter for a given HEA composition in the
    composition space of (Ti,Zr,Hf,V,Nb,Ta,Mo,W,Re,Ru) based on Hu's 2021 model (10.1016/j.actamat.2021.116800).

    Args:
        comp: A composition string which will be cast into pymatgen Composition object or ready Composition object.
        outputType: A setting to select whether the model will output a minimalistic ordered array of values (default) 
            or dictionary of labeled values. Currently implemented options are ['array', 'dict'] and the first one is 
            default.
    Returns:
        A float list representing the predicted GSF, Surf, and D parameter. Or a labeled dictionary of these output values.
    """

    assert isinstance(comp, (str, Composition)), \
        "comp must be a string or a pymatgen Composition object."
    if isinstance(comp, str):
        comp = Composition(comp)

    assert all([e.symbol in elementsSpace for e in comp.elements]), \
        "The composition must be in the composition space of (Ti,Zr,Hf,V,Nb,Ta,Mo,W,Re,Ru)."

    compList = [comp.get_atomic_fraction(e) for e in elementsSpace]
    with conversion.localconverter(default_converter):
        result = heaPredFunc(robjects.FloatVector(compList), path)
    result = list(result)
    
    assert len(result)==3

    if outputType == "array":
        return result
    elif outputType == "dict":
        return {"gfse": result[0],
                "surf": result[1],
                "dparam": result[2]}
    else:
        raise ValueError("Not recognized output type requested.")

class Model():
    def __init__(self):
        self.base = importr('base')
        self.utils = importr('utils')
        self.locfit = importr('locfit')

        with conversion.localconverter(default_converter):
            self.r = robjects.r
            self.path = str(resources.files('pqam_dparamhu2021'))
            self.r['source'](self.path+'/HEA_pred.R')
            # Initialize the models
            self.heaPredInit = robjects.globalenv['init']
            self.heaPredInit(self.path)

        # Load the prediction function
        self.heaPredFunc = robjects.globalenv['HEA_pred']

    def predict(
            self,
            comp: Union[Composition, str],
            outputType: str = "array") -> Union[dict, list]:
        """
        Predicts the GSF, Surf, and resulting D parameter for a given HEA composition in the
        composition space of (Ti,Zr,Hf,V,Nb,Ta,Mo,W,Re,Ru) based on Hu's 2021 model (10.1016/j.actamat.2021.116800).

        Args:
            comp: A composition string which will be cast into pymatgen Composition object or ready Composition object.
            outputType: A setting to select whether the model will output a minimalistic ordered array of values (default)
                or dictionary of labeled values. Currently implemented options are ['array', 'dict'] and the first one is
                default.
        Returns:
            A float list representing the predicted GSF, Surf, and D parameter. Or a labeled dictionary of these output values.
        """

        assert isinstance(comp, (str, Composition)), \
            "comp must be a string or a pymatgen Composition object."
        if isinstance(comp, str):
            comp = Composition(comp)

        assert all([e.symbol in elementsSpace for e in comp.elements]), \
            "The composition must be in the composition space of (Ti,Zr,Hf,V,Nb,Ta,Mo,W,Re,Ru)."

        compList = [comp.get_atomic_fraction(e) for e in elementsSpace]
        with conversion.localconverter(default_converter):
            result = self.heaPredFunc(robjects.FloatVector(compList), self.path)
        result = list(result)

        assert len(result)==3

        if outputType == "array":
            return result
        elif outputType == "dict":
            return {"gfse": result[0],
                    "surf": result[1],
                    "dparam": result[2]}
        else:
            raise ValueError("Not recognized output type requested.")


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

if __name__ == "__main__":
    assert len(sys.argv) == 2
    print(predict(Composition(sys.argv[1]), outputType="array"))
