import chemeos.antoine as an

# Global string variables to denote temperature units.
KELVIN = "Kelvin"
CELSIUS = "Celsius"
FARENHEIT = "Farenheit"
RANKINE = "Rankine"
# Global string variables to denote pressure units.
MMHG = "mmHg"
KPA = "kPa"
BAR = "bar"
ATM = "atm"

METHANOL = "methanol"
BENZENE = "benzene"

boublik_et_al = an.AntoineCoefficientLib(CELSIUS,MMHG,"T. Boublik, V. Fried, and E. Hala\
    , Vapour Pressures of Pure Substances, Elsevier, Amsterdam, 1973")
boublik_et_al.add_coefficients(METHANOL,8.08097,1582.271,239.726,14.9,\
83.7)

print(boublik_et_al.get_a(METHANOL,50))
print(str(boublik_et_al.get_saturation_pressure(METHANOL,50)) + " " + boublik_et_al.pressure_units)


