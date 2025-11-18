import cobra
import pandas as pd
import numpy as np
from scipy.integrate import odeint
import plotly.graph_objects as go
import pytfa
from pytfa.io import import_matlab_model, load_thermoDB
import gurobipy

# set paths to model location
path_SPOB = "C:/Users/syeghiazaryan/Desktop/metabolic models/with_thermodynamics/MAG4_restricted.xml"
path_SAOB = "C:/Users/syeghiazaryan/Desktop/metabolic models/with_thermodynamics/MAG9_restricted.xml"
path_HM = "C:/Users/syeghiazaryan/Desktop/metabolic models/with_thermodynamics/MAG1_restricted.xml"

# upload information about concentration of metabolties. Change of concentration of this metabolites will be tracked
media_community = pd.read_excel("C:/Users/syeghiazaryan/Desktop/metabolic models/with_thermodynamics/media_AA_community.xlsx",
                                sheet_name='Final')
# upload information about maximum internal flux of metabolites of each species.
media_SPOB = pd.read_excel("C:/Users/syeghiazaryan/Desktop/metabolic models/with_thermodynamics/media_AA_SPOB - Copy.xlsx", sheet_name='Final')
media_SAOB = pd.read_excel("C:/Users/syeghiazaryan/Desktop/metabolic models/with_thermodynamics/media_AA_SAOB - Copy.xlsx", sheet_name='Final')
media_HM = pd.read_excel("C:/Users/syeghiazaryan/Desktop/metabolic models/with_thermodynamics/media_AA_HM - Copy.xlsx", sheet_name='Final')
# upload database with formation ΔG of metabolites
### Load database of formation ΔG of metabolites
thermo_data = load_thermoDB("C:\\path\\to\\pytfa\\data\\thermo_data.thermodb")
thermo_data_filtered = load_thermoDB("C:\\path\\to\\pytfa\\data\\thermo_data.thermodb")

#upload models
model_SPOB = cobra.io.read_sbml_model(path_SPOB)
model_SAOB = cobra.io.read_sbml_model(path_SAOB)
model_HM = cobra.io.read_sbml_model(path_HM)

# set formation ΔG of Na and K to be 0, so in transport reacting where this ions participate ΔG will be calculated
thermo_data_filtered['metabolites']['cpd00971']['deltaGf_std'] = 0
thermo_data_filtered['metabolites']['cpd00971']['deltaGf_err'] = 0
thermo_data_filtered['metabolites']['cpd00205']['deltaGf_std'] = 0
thermo_data_filtered['metabolites']['cpd00205']['deltaGf_err'] = 0

### Filter database with thermodynamics defalut values of metabolites from entries where error on the
### metabolite's thermdynamic data didn't have value "Nil".
for metabolite in thermo_data['metabolites']:
    met_dict = thermo_data['metabolites'][metabolite]
    if met_dict['error'] != 'Nil':
        del thermo_data_filtered['metabolites'][metabolite]

### Update metabolite compartment annotation to match requirements of pyTFA tool
def annotator(cobra_model):
    for metabolite in cobra_model.metabolites:
        seed_id = metabolite.id
        metabolite.annotation['seed_id'] = seed_id[:-3]
        metabolite.charge = float(metabolite.charge)
        if metabolite.compartment == 'c0':
            metabolite.compartment = 'c'
        elif metabolite.compartment == 'c_0':
            metabolite.compartment = 'c'
        elif metabolite.compartment == 'p0':
            metabolite.compartment = 'p'
        elif metabolite.compartment == 'e0':
            metabolite.compartment = 'e'

# apply annonation addition to models
annotator(model_SPOB)
annotator(model_HM)
annotator(model_SAOB)

### Set compartments information to match requirements
new = {}
new['c'] = {'pH': 7.5, 'ionicStr': 0.25, 'symbol': 'c', 'name': 'Cytosol',
  'c_max': 0.05, 'c_min': 5e-16, 'membranePot': {'c': 0.0, 'm': 0.0, 'v': 0.0, 'x': 0.0, 'e': 150.0,
   't': 0.0, 'g': 0.0, 'r': 0.0, 'n': 0.0, 'p': 0.0}}
new['e'] = {'pH': 8.5, 'ionicStr': 0.0, 'symbol': 'e', 'name': 'Extra-organism',
  'c_max': 0.1, 'c_min': 1e-09, 'membranePot': {'c': -150.0, 'm': 0.0, 'v': 0.0, 'x': 0.0,
   'e': 0.0, 't': 0.0, 'g': 0.0, 'r': 0.0, 'n': 0.0, 'p': 0.0}}
new['p'] = {'pH': 7.5, 'ionicStr': 0.0, 'symbol': 'p', 'name': 'Periplasm', 'c_max': 0.05,
  'c_min': 1e-06, 'membranePot': {'c': 0.0, 'm': 0.0, 'v': 0.0, 'x': 0.0, 'e': 0.0, 't': 0.0,
   'g': 0.0, 'r': 0.0, 'n': 0.0, 'p': 0.0}}

# insert new compartment information into models
model_SPOB.compartments = new
model_SAOB.compartments = new
model_HM.compartments = new

# add thermodynamic constraints into models
model_SPOB_th = pytfa.ThermoModel(thermo_data_filtered, model_SPOB)
model_SAOB_th = pytfa.ThermoModel(thermo_data_filtered, model_SAOB)
model_HM_th = pytfa.ThermoModel(thermo_data_filtered, model_HM)

model_SPOB_th.prepare()
model_SAOB_th.prepare()
model_HM_th.prepare()

model_SPOB_th.convert()
model_SAOB_th.convert()
model_HM_th.convert()

# set out variables for future use
substrate_Id = list(media_community['exchange reaction ID'])
substrate_name = list(media_community['Formula'])
substrate_concentration = list(media_community['Concentration (mmol/L)'])
intrinsic_flux_SPOB = [x for x in list(media_SPOB['flux'] * -1)]
intrinsic_flux_SAOB = [x for x in list(media_SAOB['flux'] * -1)]
intrinsic_flux_HM = [x for x in list(media_HM['flux'] * -1)]

Par = {}
# define state variables: we will screen in time 3 entities: Ecoli, glucose and valine
Par['state_variable'] = {'community': substrate_name, 'SPOB': list(media_SPOB['Formula']),
                         'SAOB': list(media_SAOB['Formula']), "HM": list(media_HM['Formula'])}

# initial condition for each entity
Par['initiale_state'] = dict(zip(substrate_name, substrate_concentration))
Par['variable_model'] = ['Time'] + Par['state_variable']['community']
# Mapping to navigate between numpy arrays and state variable
Par['Id'] = {variable: index for index, variable in enumerate(Par['state_variable']['community'])}
# Choose the model
Par['model'] = {'SPOB': model_SPOB_th, 'SAOB': model_SAOB_th, 'HM': model_HM_th}

# Choose the intrinsic flux (see diapositive)
Par['intrinsic_flux'] = {'SPOB': dict(zip(media_SPOB['Formula'], intrinsic_flux_SPOB)),
                         'SAOB': dict(zip(media_SAOB['Formula'], intrinsic_flux_SAOB)),
                         "HM": dict(zip(media_HM['Formula'], intrinsic_flux_HM))}

Par['flux_id'] = {'keys': substrate_name,
                  'model': dict(zip(substrate_name, substrate_Id))}
# solution dimension
Par['Number'] = {'sol_dfba': 1201,
                 'time_point_dfba': 10}
# time step
Par['time_step_imp'] = 2

model_HM_th.solver = 'gurobi'
model_SAOB_th = 'gurobi'
model_SPOB_th = 'gurobi'

def compute_RHS(y, t, Par):
    """
    Input:
    - y: numpy array. State vector
    - t: scalar. Current time. We do not use this parameter in the function,
    but it is necessary to add it in order to function computer_RHS comply
    with scipy.integrate RHS function signature.
    - Par: dictionary. Dynamic system parameter dictionary.
    """
    output = np.zeros_like(y)
    for microbe in ['SPOB', 'SAOB', 'HM']:
        with Par['model'][microbe] as mod:
            # extract important information from state vector
            bio = y[Par['Id'][microbe]]
            # compute lower bound for substrate
            constr = {}

            for i in Par['state_variable'][microbe]:
                if y[Par['Id'][i]] < 0:
                    constr[i] = 0
                else:
                    constr[i] = max(Par['intrinsic_flux'][microbe][i],
                                    -(y[Par['Id'][i]] / (Par['time_step_imp'] * bio)))

            # affect lower bound
            for i in Par['state_variable'][microbe]:
                id_reaction = Par['flux_id']['model'][i]
                reaction = mod.reactions.get_by_id(id_reaction)
                reaction.lower_bound = constr[i]

            # compute FBA
            try:
                solution = mod.optimize()
                solut = {}

                if solution.status == 'optimal':
                    for i in Par['state_variable'][microbe]:
                        solut[i] = solution[Par['flux_id']['model'][i]]

                else:
                    for i in Par['state_variable'][microbe]:
                        solut[i] = 0
            except gurobipy.GurobiError:
                for i in Par['state_variable'][microbe]:
                    solut[i] = 0

            for i in Par['state_variable'][microbe]:
                output[Par['Id'][i]] += solut[i] * bio

    return output


y0 = np.array([float(Par['initiale_state'][k]) for k in Par['state_variable']['community']])

initial_time = 0

final_time = 2

time = np.linspace(initial_time, final_time, Par['Number']['sol_dfba'])

sol = odeint(compute_RHS, y0, time, args=(Par,))

fig = go.Figure()

# Plot variables except 'SAOB' on primary y-axis
for var in Par['state_variable']['community']:
    if var not in ['SAOB', 'SPOB', 'HM']:
        fig.add_trace(go.Scatter(
            x=time,
            y=sol[:, Par['Id'][var]],
            name=var,
            yaxis='y1',
            line=dict(width=10)
        ))

# Plot 'SAOB' on secondary y-axis
for microbe in ['SPOB', 'SAOB', 'HM']:
    fig.add_trace(go.Scatter(
        x=time,
        y=sol[:, Par['Id'][microbe]],
        name=microbe,
        yaxis='y2',
        line=dict(width=10)
    ))

# Set up layout with dual y-axes
fig.update_layout(
    width=1100,
    height=700,
    xaxis=dict(title='Time'),
    yaxis=dict(title='State Variables'),
    yaxis2=dict(title='SAOB', overlaying='y', side='right'),
    legend=dict(x=0.01, y=0.99),
    template="simple_white"
)

fig.show()


