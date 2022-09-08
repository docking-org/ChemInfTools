# Matthew Smith, 11/25/19
# Shoichet Lab: Torsion Strain Rotation Project

# In this Python script, we present the revised method for turning the
# frequency information for dihedral angles in a Torsion Library XML file into
# energy information. The purpose is to then use the new XML file, which
# preserves the original library's structure, to estimate the energy of a
# dihedral angle in a given conformation for a molecule of interest. This
# version is updated from Version 5 by the addition of a "histogram_updated"
# child element for every exact-method torsion pattern, instead of a
# "histogram_converted". The histogram_updated will just be the
# histogram_shifted with the minimum count added to any zeroes. We will
# create the histogram_converted child in the R script Add_CI.R, where we will
# also add confidence intervals to each energy estimate.

# We read in and parse the library of torsion pattern information from
# "Torsion Angle Preferences in Druglike Chemical Space: A Comprehensive
# Guide" (in the supplemental XML .txt file) by Schärfer et al. in 2013.
# We first remove the torsion patterns that contain only 3 atoms or the
# SMARTS addition "N_lp" (or &quot;N_lp&quot;). For each torsion pattern
# remaining in the library, we either directly append a "histogram_converted"
# element that converts frequencies of dihedral angles into energy (relative
# to the most frequent dihedral angle) or only convert the (value, tolerance1,
# tolerance2)-data that the authors compiled into (theta_0, beta_1,
# beta_2)-data using the math I wrote out in my notebook (assuming that around
# each peak theta_0, the energy looks like
#   1 - cosx ≈ (1/2)x^2 - (1/4!)x^4
# and then fitting coefficients to X^2 and x^4). We will dedcide which
# method to use based on the total count of observed angles (using the
# approximate method for low counts). We will read in Torsion Library
# version 2.1 as an XML file and not just parse it as a long string, so
# we will lose the comments in the original file.

# Obligatory text from TorsionLibrary 2.1:
# "This software uses the Torsion Library jointly developed by
# the University of Hamburg, Center for Bioinformatics, Hamburg,
# Germany and F. Hoffmann-La-Roche Ltd., Basel, Switzerland."

# Guba, W., Meyder, A., Rarey, M., and Hert, J. (2015). Torsion
# Library Reloaded: A New Version of Expert-Derived SMARTS Rules
# for Assessing Conformations of Small Molecules. Journal of
# Chemical Information and Modeling, DOI: acs.jcim.5b00522

# Schärfer, C., Schulz-Gasch, T., Ehrlich, H.C., Guba, W., Rarey, M.,
# Stahl, M. (2013). Torsion Angle Preferences in Drug-like Chemical
# Space: A Comprehensive Guide. Journal of Medicinal Chemistry, 56
# (6):2016-28.

import xml.etree.ElementTree as ET #For reading and writing XML files
from scipy.constants import gas_constant #Boltzmann constant in J/mol*K
from math import log #Logarithm (default is base e)
from math import pi #pi!
import numpy as np #For arrays and linear regression
from rdkit import Chem #For MolFromSmarts and GetNumAtoms


R = gas_constant / 4184
# Converting from J/mol*K to kcal/mol*K, according to Google and Wikipedia
temp = 25
# Assuming the chemical structures in the CSD came from room temperature, in C
T = temp + 273.15 #Convert from Celsius to Kelvin, based on Google/Wikipedia
N_thresh = 100
# Total count threshold for deciding whether to use the approximate method


def E_est(peak, tol_1, tol_2, n, n_2):
    # This is a function to fit the histogram of torsion angles around a
    # local peak to a Boltzmann distribution of energies. Given the local
    # peak (where the histogram is greater than 4%) and the first two
    # "tolerances" around the peak (defined as the angular distance where the
    # histogram dips below 2.5% and 1.5% on each side of the peak), we estimate
    # the energies at the different tolerances relative to the peak. We then
    # fit these five points to a quadratic function and a quartic function
    # of the dihedral angle using linear regression and report the location
    # of the minimum (theta_0 = peak) and coefficients from the regression.
    # We arrive at using x^2 and x^4 because
    #   1 - cosx ≈ (1/2)x^2 - (1/4!)x^4
    # is a much better approximation than
    #   1 - cosx ≈ (1/2)x^2
    # when the angles are large (at least up to about 90 degrees). We really
    # would be fitting to
    #   E(x) = A(1 - cos(omega*x)),
    # which has 2 parameters (A, omega), so we fit different parameters
    # to x^2 and x^4 (absorbing the signs and coefficients from the Taylor
    # Expansion).

    # We use n for the number of peaks for the given bond and n_2 for
    # the number of second tolerances that differ from first
    # tolerances. We use these in our estimates of the peak probability,
    # which we then use to find the energies of tol_1 and tol_2.

    theta_0 = peak #The most probable angle
    # Since we are subtracting the theta_0 from the other angles, we don't
    # need to calculate them explicitly as theta_0 +- tol_1 or tol_2

    P_tot = 0.04*n + 0.05*n_2 + 0.03*n #Total probability mass. See below
    if P_tot > 1:
        # Here, the probability masses contributed by the all the peaks add up
        # to > 1, which can happen if the Torsion Library authors do not
        # strictly use P(Peak) >= 0.04, P(tol_1) = 0.025, and P(tol_2) = 0.015.
        # We then rescale these probabilities by the total mass and do not
        # attempt to make a better estimate for the peak probability
        P_0 = 0.04 / P_tot
        P_1 = 0.025 / P_tot
        P_2 = 0.015 / P_tot
    else:
        # Here we can directly estimate the peak probability. We assume that
        # each peak has the same probability, but we need to account for the
        # probability mass taken up by the points at the different tolerance
        # levels. Each peak has 2 points that define its second tolerance,
        # and each one has a probability of 0.015. Wherever there is a first
        # tolerance differing from the second tolerance, we get another two
        # points with probability mass, each one now with probability 0.025.
        # We have n as the number of peaks for a given torsion pattern and
        # n_2 as the number of first tolerances that differ from the second
        # tolerances for the same peak for a given torsion pattern. Then the
        # probability for each peak is, with a total probability of 1, is:
        #   P(peak) = (1 - 2*0.015*n - 2*0.025*n_2) / n
        #           = 1/n - 0.03 - 0.05*(n_2/n)
        P_0 = 1/n - 0.03 - 0.05 * (n_2/n)
        P_1 = 0.025 #This and P_2 do not need to be rescaled
        P_2 = 0.015

    # Invert the Boltzmann weights and set E_0 = 0. See notes
    E_1 = -R*T*log(P_1) + R*T*log(P_0) #In TEU
    E_2 = -R*T*log(P_2) + R*T*log(P_0) #In TEU

    # Fit the slopes, in degrees, not radians!!!
    if tol_1 == tol_2:
        # If this is the case, then we only have two points to fit two
        # slope coefficients. We will use the new method of using
        #   E(x) = A * [1-cos(omega*x)]
        #        ≈ A * [(1/2)*(omega*x)^2 - (1/24)*(omega*x)^4]
        #        = (A*omega^2)/2 * x^2 - (A*omega^4)/24 * x^4
        #        = beta_1 * x^2 + beta_2 * x^4
        # and making tol_1 = tol_2 = tol be at the maximum value. We
        # need the derivative of E:
        #   E'(x) = beta_1 * 2x + beta_2 * 4x^3
        #         = 2x * (beta_1 + 2*beta_2*x^2)
        # The maximum value (tol, E_2) satisfies
        #   0 = E'(tol) = 2*tol * (beta_1 + 2*beta_2*tol^2)
        #   => beta_1 = -2*beta_2*tol^2
        # and
        #   E_2 = E(tol) = beta_1 * tol^2 + beta_2 * tol^4
        #       = -2*beta_2*tol^2 * tol^2 + beta_2 * tol^4
        #       = -2 * beta_2 * tol^4 + beta_2 * tol^4
        #       = -beta_2 * tol^4
        # So
        #   beta_2 = - E_2 / tol^4
        #   beta_1 = 2 * E_2 / tol^2
        # With E_2 > 0, this means beta_1 > 0 and bet_2 < 0, as expected.
        # This works when x is in degrees or radians.

        beta_1 = 2*E_2 / tol_2**2 #In TEU/deg^2
        beta_2 = -E_2 / tol_2**4 #In TEU/deg^4
    else:
        # When we have different first and second tolerances, then we can
        # just use regression to fit beta_1 and beta_2, now in degrees!
        # The regression coefficients will absorb the angle-unit conversions
        x_1 = np.array([(-tol_2)**2, (-tol_1)**2, 0, tol_1**2, tol_2**2])
        x_2 = np.array([(-tol_2)**4, (-tol_1)**4, 0, tol_1**4, tol_2**4])
        y = np.array([E_2, E_1, 0, E_1, E_2])
        X_mat = np.vstack([np.ones(len(x_1)), x_1, x_2]).T #Design matrix
        beta_0, beta_1, beta_2 = np.linalg.lstsq(X_mat, y, rcond=None)[0]
        # We won't need beta_0

    return theta_0, beta_1, beta_2


# Import the XML file, using this as a guide:
# https://docs.python.org/3/library/xml.etree.elementtree.html
tree = ET.parse('TorLibv21WCSDBins.xml')
root = tree.getroot()


# We can loop directly over each of the (currently) 514 torsion rules!
# Thank you ET!
# We use the iter() function to iterate over all of the torsionRules,
# regardless of whether they are children of the hierarchy class (like for
# most of GG), grandchildren of the hierarchy class (like for the "aro aro"
# subclass of GG), or great-grandchildren of the hierarchy class (like for
# all of the specific classes)
# To remove an element, we need to call the .remove method on its parent.
# ET does not keep track of each element's parent by default, so we need
# to create a dictionary matching each child element to its parent, based on
#https://stackoverflow.com/questions/2170610/access-elementtree-node-parent-node
parent_map = {c:p for p in root.iter() for c in p}

# First we determine whether the SMARTS pattern for the torsion pattern
# contains the SMARTS addition "N_lp" (or &quot;N_lp&quot;) or only 3
# atoms. We remove these problematic torsion patterns. We need to loop over
# the library multiple times, because the .remove method sometimes does not
# remove multiple child nodes from the same
found_TP = True #Initialize whether the loop found a bad torsion pattern
while(found_TP):
    found_TP = False #Switch to false before find any
    for torsionRule in root.iter("torsionRule"):
        smarts = torsionRule.get("smarts")
        if "N_lp" in smarts: #Nifty Python string-checking!
            found_TP = True #Now that the loop found a bad torsion pattern
            # print(parent_map[torsionRule]) #Debugging
            parent_map[torsionRule].remove(torsionRule)
            # Remove from the parent element
            # print(smarts) #Debugging
        else:
            # Cannot use GetNumAtoms if "N_lp" is present, because then
            # the SMARTS pattern will not work with Chem
            atoms = Chem.MolFromSmarts(smarts).GetNumAtoms() #Number of atoms
            if atoms < 4:
                found_TP = True #Now that the loop found a bad torsion pattern
                # print(parent_map[torsionRule]) #Debugging
                parent_map[torsionRule].remove(torsionRule)
                # Remove from the parent element
                # print(smarts + ", " + atoms) #Debugging

for torsionRule in root.iter("torsionRule"): #Loop again
    # Next we determine whether the total count of the angles is less than
    # the threshold
    N = 0 #Initialize
    for bin in torsionRule.find("histogram_shifted").findall("bin"):
        # Find the histogram_shifted for each torsion rule, and then
        # every bin within
        N += int(bin.get("count")) #Add the count from each bin

    if N < N_thresh: #Use the approximate method
        # First calculate n and n_2 for E_est
        n = 0 #Initialize total number of peaks
        n_2 = 0 #Initialize the number of peaks where tol_1 != tol_2
        for angle in torsionRule.find("angleList").findall("angle"):
            # Find every angle in the list of peaks and tolerances
            n += 1
            if float(angle.get("tolerance1")) != float(angle.get("tolerance2")):
                # If the first and second tolerances are different
                n_2 += 1

        # Now that we have values for n and n_2 we can use E_est
        for angle in torsionRule.find("angleList").findall("angle"):
            p = float(angle.get("value")) #The location of the peak
            t_1 = float(angle.get("tolerance1")) #First tolerance
            t_2 = float(angle.get("tolerance2")) #Second tolerance
            theta_0, beta_1, beta_2 = E_est(p, t_1, t_2, n, n_2)
            angle.set("theta_0", str(theta_0)) #Set theta_0
            angle.set("beta_1", str(beta_1)) #Set beta_1
            angle.set("beta_2", str(beta_2)) #Set beta_2

        torsionRule.set("method", "approximate")
        # Add an attribute saying we used the approximate method

    else: #If N >= N_thresh, use the exact method with histogram_shifted
        # First we make a list with the values of the counts
        counts = [] #Initialize
        for bin in torsionRule.find("histogram_shifted").findall("bin"):
            c = int(bin.get("count"))
            counts.append(c)
        # Next we need to set all the zero counts to the minimum count
        m = min(i for i in counts if i > 0)
        # The minimum non-zero count, using a generator expression
        for i in range(len(counts)):
            if counts[i] == 0:
                counts[i] = m #Change the zero count
        # Now we create the histogram_updated object. We will calculate
        # the energy estimates and confidence intervals in the Add_CI.R
        # R script
        h = ET.SubElement(torsionRule, "histogram_updated") #Initialize
        for i in range(len(counts)):
            bin = ET.SubElement(h, "bin") #Create a new bin subelement
            bin.set("count", str(counts[i])) #Set the count attribute
        torsionRule.set("method", "exact")
        # Add an attribute saying we used the exact method


# Now we write the completed XML file
tree.write("TL_2.1_pre.xml")
# We write the suffix "_pre" because we need to add the preamble by hand
