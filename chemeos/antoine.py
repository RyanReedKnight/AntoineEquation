import numpy as np

_LOG10 = "log base 10"
_LN = "natural log"

# Global variables used as keys used in AntoineCoefficientLib
_COEFF_A_KEY = 'A'
_COEFF_B_KEY = 'B'
_COEFF_C_KEY = 'C'

# Temperature limit keys
_UPPER_TEMPERATURE_LIMIT_KEY = "upper temperature limit"
_LOWER_TEMPERATURE_LIMIT_KEY = "lower temperature limit"

# Global string variables to denote pressure units.

class SpeciesNotFound(ValueError):
    """Error raised when search for species_name in dictionary fails."""

class TemperatureOutOfRange(ValueError):
    """Error raised when temperature argument is out of range."""

def saturation_pressure(coeff_a_val,coeff_b_val,coeff_c_val, temp,is_log10 = True):
    """Given Antoine values A, B and C, and the temperature; returns the saturation pressure.

    Parameters
        - coeff_a_val is Antoine coefficient A.
        - coeff_b_val is Antoine coefficient B.
        - coeff_c_val is Antoine coe
        - temp is the temperature
    Note: It is the users responsibility to ensure units are consistent.
    """

    exponent = coeff_a_val - coeff_b_val/(coeff_c_val + temp)

    if is_log10 is True:

        return 10**exponent

    return np.e**exponent

def _check_temperature(temperature,dct):
        """Checks temperature to ensure it is within range.

        Parameters
            - temperature
            - dct: looks for the temperature range in the dictionary
                and checks if temperature is within the range. Returns
                True bool is it is in the range, False if it is not.
        """

        return bool(temperature >= dct[_LOWER_TEMPERATURE_LIMIT_KEY] and\
            temperature <= dct[_UPPER_TEMPERATURE_LIMIT_KEY])


class AntoineCoefficientLib:
    """ A class to store a dictionary of Antoine coefficients associated with chemical\
    species and perform relavent calculations.

    Fields
        source : str
            Stores a string which contains a citation for the source of the coefficients,
            instances of the class should be associated with one source.
        temperature_units: str
            The temperature units used for the coefficients, it is the users responsibility
            to keep units consistent.
        pressure_units : str
            The pressure units used for the coefficients.
        antoine_coeff_lib : dictionary
            {species_name:{'A':coeff_a_val,
            'B': coeff_b_val,
            'C':coeff_c_val,
            "upper temperature limit":upper_temperature_limit_val,
            "lower temperature limit":lower_temperature_limit_val}}
        is_log10 : boolean
            If true calculations done with log base 10, otherwise natural log is used.

    """
    def __init__(self,temperature_units,pressure_units,source,is_log10 = True):

        self.source = source
        self.temperature_units = temperature_units
        self.pressure_units = pressure_units
        self.antoine_coeff_lib = {}
        self.is_log10 = is_log10

    def _log(self,val):
        """ If _log10 is true, returns log base 10 of value, otherwise returns\
         the natural log of value."""

        if self.is_log10 is True:

            return np.log10(val)

        return np.ln(val)

    def add_coefficients(self,species_name,coeff_a_val,coeff_b_val,coeff_c_val,\
        lower_temperature_limit_val,upper_temperature_limit_val):
        """Adds a set of coefficients and parameters associeted with a chemical species.

        Parameters

            - species_name: A string with the name of the chemical species,
                S/B named the same way as the source.
            - coeff_a_val: A float containing the value of Antoine coefficient A
                for the chemical species.
            -coeff_b_val: A float containing the value of Antoine coefficient B
                for the chemical species.
            -coeff_c_value: A float containing the value of Antoine coefficient C
                for the chemical species.
            -lower_temperature_limit_val: A float containing the lower temperature limit
                of the parameters
            -upper_temperature_limit_val: A float containing the upper temperature limit
                of the parameters.

        Notes

            - The coefficients added S/B from the source associated with the instance of
                the class.
            - The user must ensure that the units are the same as those associated with the class,
                this ought not to be an issue if the user is keeping the source consistent.

        """

        new_item = {_COEFF_A_KEY:coeff_a_val,_COEFF_B_KEY:coeff_b_val,\
        _COEFF_C_KEY:coeff_c_val,_LOWER_TEMPERATURE_LIMIT_KEY\
        :lower_temperature_limit_val,_UPPER_TEMPERATURE_LIMIT_KEY\
        :upper_temperature_limit_val}

        if species_name in self.antoine_coeff_lib:

            self.antoine_coeff_lib[species_name].append(new_item)

        else:

            self.antoine_coeff_lib[species_name] = [new_item]


    def get_coeffs(self,species_name,temperature):
        """Returns a dictionary with antoine cofficients A, B and C; lower \
        tempurature limit, and upper temperature limit.

        Dictionary keys

            - key to get coefficient A is 'A'.
            - key to get coefficient A is 'B'.
            - key to get coefficient A is 'C'.
            - key to get the lower temperature limit is "lower temperature limit"
            - key to get the upper temperature limit is "upper temperature limit"

        """

        if species_name not in self.antoine_coeff_lib:

            raise SpeciesNotFound(species_name + " is not in this dictionary")

        coeffs = self.antoine_coeff_lib[species_name]

        for dct in coeffs:

            if _check_temperature(temperature,dct) is True:

                return dct

        raise TemperatureOutOfRange("There is no availible set of coefficients\
         with a temperature range " + str(temperature) + " " + str(self.temperature_units) +\
         " is within the bounds of.")

    def get_a(self,species_name,temperature):
        """ Returns coefficnet A if the temperature is not out of range."""

        return self.get_coeffs(species_name,temperature)[_COEFF_A_KEY]

    def get_b(self,species_name,temperature):
        """ Returns coefficnet B is the temperature is not out of range."""

        return self.get_coeffs(species_name,temperature)[_COEFF_B_KEY]

    def get_c(self,species_name,temperature):
        """ Returns coefficnet C is the temperature is not out of range."""

        return self.get_coeffs(species_name,temperature)[_COEFF_C_KEY]

    def get_saturation_pressure(self,species_name,temperature):
        """If species is in dictionary and temperature is in range,\
         returns saturation pressure."""

        if species_name not in self.antoine_coeff_lib:
            raise SpeciesNotFound()

        coeffs = self.get_coeffs(species_name,temperature)

        if _check_temperature(temperature,coeffs) is False:

            raise TemperatureOutOfRange()

        return saturation_pressure(coeffs[_COEFF_A_KEY],coeffs[_COEFF_B_KEY],\
            coeffs[_COEFF_C_KEY],temperature,self.is_log10)


if __name__ == "__main__":
    print("antoine is main.")
    