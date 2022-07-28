
class ModelParameterConfig:
  def __init__(self, save_every_n_steps:int, n:int, time_step:float, steps:int, 
    ecm_k_elast:float, ecm_d_dumping:float, ecm_mass:float, ecm_gel_concentration:float,
    boundary_coords:list, boundary_disp_rates:list,boundary_disp_rates_parallel:list,
    poisson_dirs:list,
    allow_boundary_elastic_movement:list, 
    boundary_stiffness:list,boundary_dumping:list,
    clamp_agent_touching_boundary:list,allow_agent_sliding:list,
    ecm_ecm_equilibrium_distance:float,ecm_boundary_interaction_radius:float,
    ecm_boundary_equilibrium_distance:float,
    include_fiber_alignment:int,ecm_orientation_rate:float,
    oscillatory_shear_assay:bool,oscillatory_amplitude:float,oscillatory_freq:float,oscillatory_w:float,
    include_diffusion:bool,n_species:int,diffusion_coeff_multi:float,
    boundary_conc_init_multi:list,boundary_conc_fixed_multi:list,
    init_ecm_concentration_vals:list,
    include_vascularization:bool, init_vascularization_concentration_vals:list):
    
    self.SAVE_EVERY_N_STEPS = save_every_n_steps      
    self.N = n
    self.TIME_STEP = time_step
    self.STEPS = steps
    self.ECM_K_ELAST = ecm_k_elast            
    self.ECM_D_DUMPING = ecm_d_dumping     
    self.ECM_MASS = ecm_mass
    self.ECM_GEL_CONCENTRATION = ecm_gel_concentration 
    self.BOUNDARY_COORDS = boundary_coords
    self.BOUNDARY_DISP_RATES = boundary_disp_rates
    self.BOUNDARY_DISP_RATES_PARALLEL = boundary_disp_rates_parallel
    self.POISSON_DIRS = poisson_dirs
    self.ALLOW_BOUNDARY_ELASTIC_MOVEMENT = allow_boundary_elastic_movement
    self.BOUNDARY_STIFFNESS = boundary_stiffness
    self.BOUNDARY_DUMPING = boundary_dumping
    self.CLAMP_AGENT_TOUCHING_BOUNDARY = clamp_agent_touching_boundary
    self.ALLOW_AGENT_SLIDING = allow_agent_sliding
    self.ECM_ECM_EQUILIBRIUM_DISTANCE = ecm_ecm_equilibrium_distance
    self.ECM_BOUNDARY_INTERACTION_RADIUS = ecm_boundary_interaction_radius
    self.ECM_BOUNDARY_EQUILIBRIUM_DISTANCE = ecm_boundary_equilibrium_distance
    self.INCLUDE_FIBER_ALIGNMENT = include_fiber_alignment
    self.ECM_ORIENTATION_RATE = ecm_orientation_rate
    self.OSCILLATORY_SHEAR_ASSAY = oscillatory_shear_assay
    self.OSCILLATORY_AMPLITUDE = oscillatory_amplitude
    self.OSCILLATORY_FREQ = oscillatory_freq
    self.OSCILLATORY_W = oscillatory_w 
    self.INCLUDE_DIFFUSION = include_diffusion
    self.N_SPECIES = n_species 
    self.DIFFUSION_COEFF_MULTI = diffusion_coeff_multi                               
    self.BOUNDARY_CONC_INIT_MULTI = boundary_conc_init_multi                            
    self.BOUNDARY_CONC_FIXED_MULTI = boundary_conc_fixed_multi                             
    self.INIT_ECM_CONCENTRATION_VALS = init_ecm_concentration_vals
    self.INCLUDE_VASCULARIZATION = include_vascularization    
    self.INIT_VASCULARIZATION_CONCENTRATION_VALS = init_vascularization_concentration_vals



