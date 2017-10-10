Constants <- list(
    MAX_SIMULATION_SPEED = 1000,
    # Measures (in ms) the amount of time considered to respond against an incident
    INCIDENT_ACTION_THRESHOLD = 10000,    
    # The minimum time that a simulation is considered to last to not being considered as an outlier
    OUTLIER_REAL_TIME_THRESHOLD = 20000,
    INTERACTIONS_THRESHOLD = 10,
    MIN_TIME_BETWEEN_INTERACTIONS = 0.5,
    MAX_TIME_BETWEEN_INTERACTIONS_THRESHOLD = 8, # seconds
    MIN_INTERACTIONS_PER_SIMULATION = 10,    
    INTERACTION_NAMES = c("SELECT\nUAV",
                          "SET\nUAV\nSPEED",
                          "SET\nSIMULATION\nSPEED",
                          "CHANGE\nUAV\nPATH",
                          "SET\nWAYPOINT\nTABLE",
                          "SET\nCONTROL\nMODE")
  )