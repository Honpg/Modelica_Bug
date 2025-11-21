package HVAC_FMU "AHU Complete – RL-READY with External Weather Data"
  import Modelica;
  import Buildings;

  // ========= Medium definitions =========
  package MediumA = Buildings.Media.Air(
    extraPropertiesNames = {"CO2"},
    C_nominal = {1.519e-3});
  package MediumW = Buildings.Media.Water;

  /////////////////////////////////////////////////////////////////////////////
  // 1) Chiller with manual control - ENERGY BALANCED
  /////////////////////////////////////////////////////////////////////////////

  model Carnot_TEva_NoWarn
    "Chiller with TEva setpoint and manual on/off control - Energy balanced"
    extends Buildings.Fluid.Chillers.BaseClasses.PartialCarnot_T(
      redeclare package Medium1 = Buildings.Media.Water,
      redeclare package Medium2 = Buildings.Media.Water,
      final COP_is_for_cooling = true,
      final QCon_flow_nominal = -QEva_flow_nominal*(1 + COP_nominal)/COP_nominal,
      PEle(y = -QEva_flow/COP),
      redeclare Buildings.Fluid.HeatExchangers.HeaterCooler_u con(
        final from_dp = from_dp1,
        final dp_nominal = dp1_nominal,
        final linearizeFlowResistance = linearizeFlowResistance1,
        final deltaM = deltaM1,
        final tau = tau1,
        final T_start = T1_start,
        final energyDynamics = energyDynamics,
        final homotopyInitialization = homotopyInitialization,
        final Q_flow_nominal = QCon_flow_nominal),
      redeclare Buildings.Fluid.HeatExchangers.SensibleCooler_T eva(
        final from_dp = from_dp2,
        final dp_nominal = dp2_nominal,
        final linearizeFlowResistance = linearizeFlowResistance2,
        final deltaM = deltaM2,
        final QMin_flow = QEva_flow_min,
        final tau = tau2,
        final T_start = T2_start,
        final energyDynamics = energyDynamics,
        final homotopyInitialization = homotopyInitialization));

    parameter Modelica.Units.SI.HeatFlowRate QEva_flow_min(max = 0) = -Modelica.Constants.inf
      "Maximum cooling capacity (negative)";
    parameter Modelica.Units.SI.Temperature TEvaSet_cold = 279.15
      "Evaporator setpoint when chiller at full capacity (cold) [K]";
    parameter Modelica.Units.SI.Temperature TEvaSet_off = 298.15
      "Evaporator setpoint when chiller is off (warm) [K]";

    Modelica.Blocks.Interfaces.RealInput TSet(unit = "K")
      "Base evaporator leaving water temperature setpoint"
      annotation(Placement(transformation(extent = {{-140, 70}, {-100, 110}})));
    Modelica.Blocks.Interfaces.RealInput u(min = 0, max = 1)
      "Chiller capacity control signal [0-1]"
      annotation(Placement(transformation(extent = {{-140, 40}, {-100, 80}})));

  protected
    Modelica.Blocks.Math.Add QCon_flow_internal(final k1 = -1)
      "Heat added to condenser"
      annotation(Placement(transformation(extent = {{-80, 30}, {-60, 50}})));
    Modelica.Blocks.Sources.RealExpression yCon(
      y = QCon_flow_internal.y/QCon_flow_nominal)
      "Normalized condenser heat flow rate"
      annotation(Placement(transformation(extent = {{-40, 30}, {-20, 50}})));
    Modelica.Blocks.Sources.RealExpression TEvaSet_eff(
      y = u*TEvaSet_cold + (1 - u)*TEvaSet_off)
      "Effective evaporator setpoint modulated by control signal"
      annotation(Placement(transformation(extent = {{-80, 70}, {-60, 90}})));

  initial equation
    assert(QEva_flow_nominal < 0, "Parameter QEva_flow_nominal must be negative.");

  equation
    connect(TEvaSet_eff.y, eva.TSet) annotation(
      Line(points = {{-59, 80}, {-50, 80}, {-50, 90}, {28, 90}, {28, -52}, {12, -52}},
           color = {0, 0, 127}));
    connect(eva.Q_flow, QEva_flow) annotation(
      Line(points = {{-11, -52}, {-40, -52}, {-40, -90}, {110, -90}},
           color = {0, 0, 127}));
    connect(QCon_flow_internal.y, QCon_flow) annotation(
      Line(points = {{-59, 40}, {-52, 40}, {-52, 80}, {80, 80}, {80, 90}, {110, 90}},
           color = {0, 0, 127}));
    connect(QCon_flow_internal.u1, eva.Q_flow) annotation(
      Line(points = {{-82, 46}, {-90, 46}, {-90, -52}, {-11, -52}},
           color = {0, 0, 127}));
    connect(QCon_flow_internal.u2, PEle.y) annotation(
      Line(points = {{-82, 34}, {-88, 34}, {-88, 20}, {72, 20}, {72, 0}, {61, 0}},
           color = {0, 0, 127}));
    connect(yCon.y, con.u) annotation(
      Line(points = {{-19, 40}, {-16, 40}, {-16, 66}, {-12, 66}},
           color = {0, 0, 127}));
    annotation(
      Icon(
        coordinateSystem(
          preserveAspectRatio = false,
          extent = {{-100, -100}, {100, 100}}),
        graphics = {
          Text(
            extent = {{-148, 156}, {-92, 114}},
            textColor = {0, 0, 127},
            textString = "TEva")}));
  end Carnot_TEva_NoWarn;

  ///////////////////////////////////////////////
  // 2) AHU model - RL-READY with External Weather Data
  ///////////////////////////////////////////////

  model AHU_FMU_Core
    "Complete AHU - RL-READY with External Weather Data Inputs"
    extends Modelica.Icons.Example;

    // ============ PARAMETERS ============
    parameter Modelica.Units.SI.MassFlowRate mAir_flow_nominal = 0.24;
    parameter Modelica.Units.SI.MassFlowRate mWat_flow_nominal = 0.25;
    parameter Modelica.Units.SI.Volume Vzone = 168;
    parameter Real COP_chiller_nominal = 3.0;
    parameter Modelica.Units.SI.Power Pfan_nominal = 70;
    parameter Modelica.Units.SI.Power PfanEA_nominal = 50;
    parameter Modelica.Units.SI.Power PpumpCW_nominal = 180;
    parameter Modelica.Units.SI.Power Q_internal = 1000;
    parameter Modelica.Units.SI.HeatFlowRate QEva_flow_nominal = -3500;
    parameter Modelica.Units.SI.Power P_chiller_nominal = abs(QEva_flow_nominal)/COP_chiller_nominal;
    parameter Modelica.Units.SI.ThermalConductance UA_nominal = 1000;
    parameter Real r_nominal = 2.0/3.0;
    parameter Integer nEle = 6;
    parameter Integer nOccMax = 10;
    parameter Real moisture_gen_per_person = 50.0/3600.0/1000.0;
    parameter Real CO2_gen_per_person = 0.0104/3600;
    parameter Real CO2_outdoor_ppm = 400 "Outdoor CO2 concentration [ppm]";
    parameter Real CO2_target_full_ppm = 1000;
    parameter Modelica.Units.SI.Temperature T_coil_dehum = 279.15;
    parameter Modelica.Units.SI.Temperature TEva_cold = 279.15;
    parameter Modelica.Units.SI.Temperature TEva_off = 298.15;

    // ============ SYSTEM ============
    inner Modelica.Fluid.System system(
      energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial,
      massDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial,
      allowFlowReversal = true,
      p_ambient = 101325,
      T_ambient = 301.15,
      m_flow_small = 1e-6,
      dp_small = 10) annotation(
        Placement(transformation(origin = {-356, 18},
                                 extent = {{-10, -10}, {10, 10}})));

    // ============ WEATHER INPUTS ============
    Modelica.Blocks.Interfaces.RealInput TDryBul_in(start = 301.15)
      "Dry bulb temperature [K]"
      annotation(Placement(transformation(origin = {-366, 122},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {-320, 114},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealInput relHum_in(start = 0.6)
      "Relative humidity [0-1]"
      annotation(Placement(transformation(origin = {-366, 108},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {-320, 100},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealInput pAtm_in(start = 101325)
      "Atmospheric pressure [Pa]"
      annotation(Placement(transformation(origin = {-366, 136},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {-320, 128},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealInput HGloHor_in(start = 0)
      "Global horizontal radiation [W/m2]"
      annotation(Placement(transformation(origin = {-366, 58},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {-320, 50},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealInput HDifHor_in(start = 0)
      "Diffuse horizontal radiation [W/m2]"
      annotation(Placement(transformation(origin = {-366, 72},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {-320, 64},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealInput winSpe_in(start = 2.0)
      "Wind speed [m/s]"
      annotation(Placement(transformation(origin = {-366, 86},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {-320, 78},
                                              extent = {{-10, -10}, {10, 10}})));

    // ============ CONTROL INPUTS - RL ACTIONS ============
    Modelica.Blocks.Interfaces.RealInput uFan(start = 1)
      "Fan speed control [0-1]"
      annotation(Placement(transformation(origin = {-286, 168},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {-260, 156},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealInput uOA(start = 1)
      "Outside air damper [0-1]"
      annotation(Placement(transformation(origin = {-346, 168},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {-320, 156},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealInput uChiller(start = 1)
      "Chiller capacity [0-1]"
      annotation(Placement(transformation(origin = {-326, 168},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {-300, 156},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealInput uHeater(start = 1)
      "Heater capacity [0-1] - RL Control"
      annotation(Placement(transformation(origin = {-366, 168},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {-340, 156},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealInput occupancy(start = 1)
      "Occupancy fraction [0-1]"
      annotation(Placement(transformation(origin = {-266, 168},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {-240, 156},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealInput uFanEA(start = 1)
      "Exhaust fan speed [0-1]"
      annotation(Placement(transformation(origin = {-306, 168},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {-280, 156},
                                              extent = {{-10, -10}, {10, 10}})));

    // ============ OUTPUTS - RL OBSERVATIONS ============
    Modelica.Blocks.Interfaces.RealOutput T_SA
      "Supply air temperature [K]"
      annotation(Placement(transformation(origin = {318, 168},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {330, 140},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealOutput RH_SA
      "Supply air relative humidity [0-1]"
      annotation(Placement(transformation(origin = {318, 148},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {330, 120},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealOutput Vdot_SA
      "Supply air volume flow rate [m3/s]"
      annotation(Placement(transformation(origin = {318, 128},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {330, 100},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealOutput T_zone
      "Zone temperature [K]"
      annotation(Placement(transformation(origin = {318, 108},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {330, 80},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealOutput RH_zone
      "Zone relative humidity [0-1]"
      annotation(Placement(transformation(origin = {318, 88},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {330, 60},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealOutput CO2_zone_ppm
      "Zone CO2 concentration [ppm]"
      annotation(Placement(transformation(origin = {318, 68},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {330, 40},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealOutput P_fan
      "Supply fan power [W]"
      annotation(Placement(transformation(origin = {278, 168},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {310, 160},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealOutput Q_chiller
      "Chiller cooling rate [W]"
      annotation(Placement(transformation(origin = {298, 168},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {330, 160},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealOutput P_chiller
      "Chiller power [W]"
      annotation(Placement(transformation(origin = {218, 168},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {248, 160},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealOutput P_fanEA
      "Exhaust fan power [W]"
      annotation(Placement(transformation(origin = {318, 48},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {330, 20},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealOutput P_pump
      "Chilled water pump power [W]"
      annotation(Placement(transformation(origin = {238, 168},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {270, 160},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealOutput P_total
      "Total system power [W]"
      annotation(Placement(transformation(origin = {258, 168},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {290, 160},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealOutput Q_heater
      "Heater heat output [W]"
      annotation(Placement(transformation(origin = {318, 28},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {330, 0},
                                              extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Interfaces.RealOutput T_SA_afterCooling
      "SA temperature after cooling coil [K]"
      annotation(Placement(transformation(origin = {318, 8},
                                          extent = {{-10, -10}, {10, 10}}),
                           iconTransformation(origin = {330, -20},
                                              extent = {{-10, -10}, {10, 10}})));

    // ============ COMPONENTS ============
    Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(
      filNam = Modelica.Utilities.Files.loadResource(
        "modelica://Buildings/Resources/weatherdata/USA_CA_San.Francisco.Intl.AP.724940_TMY3.mos"),
      TDryBulSou = Buildings.BoundaryConditions.Types.DataSource.Input,
      relHumSou  = Buildings.BoundaryConditions.Types.DataSource.Input,
      pAtmSou    = Buildings.BoundaryConditions.Types.DataSource.Input,
      winSpeSou  = Buildings.BoundaryConditions.Types.DataSource.Input,
      HSou       = Buildings.BoundaryConditions.Types.RadiationDataSource.Input_HGloHor_HDifHor)
      annotation(Placement(transformation(origin = {-304, 96},
                                          extent = {{-28, -28}, {28, 28}})));

    Buildings.Fluid.Sources.Outside ambAir(
      redeclare package Medium = MediumA,
      nPorts = 2,
      use_C_in = false,
      C = {CO2_outdoor_ppm*1e-6})
      "Outdoor air boundary with fixed CO2 concentration"
      annotation(Placement(transformation(origin = {-248, 96},
                                          extent = {{-12, -12}, {12, 12}})));

    Modelica.Blocks.Continuous.FirstOrder filtOA(
      T = 90,
      initType = Modelica.Blocks.Types.Init.SteadyState)
      annotation(Placement(transformation(origin = {-228, 133},
                                          extent = {{-6, -6}, {6, 6}})));

    Modelica.Blocks.Nonlinear.Limiter limOA(
      uMin = 0.0,
      uMax = 1.0)
      annotation(Placement(transformation(origin = {-208, 133},
                                          extent = {{-6, -6}, {6, 6}})));

    Modelica.Blocks.Continuous.FirstOrder filtRA(
      T = 90,
      initType = Modelica.Blocks.Types.Init.InitialState,
      y_start = 0.7)
      annotation(Placement(transformation(origin = {162, -20},
                                          extent = {{-6, -6}, {6, 6}})));

    Modelica.Blocks.Continuous.FirstOrder filtEA(
      T = 90,
      initType = Modelica.Blocks.Types.Init.SteadyState)
      annotation(Placement(transformation(origin = {198, -73},
                                          extent = {{-6, -6}, {6, 6}})));

    Modelica.Blocks.Nonlinear.Limiter limEA(
      uMin = 0.0,
      uMax = 1.0)
      annotation(Placement(transformation(origin = {218, -73},
                                          extent = {{-6, -6}, {6, 6}})));

    Modelica.Blocks.Continuous.FirstOrder filtValve(
      T = 120,
      initType = Modelica.Blocks.Types.Init.InitialState,
      y_start = 0.5)
      annotation(Placement(transformation(origin = {-38, 62},
                                          extent = {{-6, -6}, {6, 6}})));

    Modelica.Blocks.Continuous.FirstOrder filtHeater(
      T = 90,
      initType = Modelica.Blocks.Types.Init.InitialState,
      y_start = 0.0)
      annotation(Placement(transformation(origin = {-56, 136},
                                          extent = {{-6, -6}, {6, 6}})));

    Buildings.Fluid.Actuators.Dampers.Exponential damOA(
      redeclare package Medium = MediumA,
      m_flow_nominal = mAir_flow_nominal,
      dpDamper_nominal = 60,
      dpFixed_nominal = 10,
      allowFlowReversal = false,
      use_constant_density = true,
      from_dp = false,
      linearized = true)
      annotation(Placement(transformation(origin = {-216, 96},
                                          extent = {{-10, -10}, {10, 10}})));

    Buildings.Fluid.Sensors.MassFlowRate senOA(
      redeclare package Medium = MediumA,
      allowFlowReversal = false)
      annotation(Placement(transformation(origin = {-186, 96},
                                          extent = {{-10, -10}, {10, 10}})));

    Buildings.Fluid.FixedResistances.PressureDrop preFil(
      redeclare package Medium = MediumA,
      m_flow_nominal = mAir_flow_nominal,
      dp_nominal = 30,
      allowFlowReversal = false,
      from_dp = false,
      linearized = true)
      annotation(Placement(transformation(origin = {-156, 96},
                                          extent = {{-10, -10}, {10, 10}})));

    // ===== FIXED: mixT portFlowDirection_* = Bidirectional =====
    Buildings.Fluid.FixedResistances.Junction mixT(
      redeclare package Medium = MediumA,
      m_flow_nominal = mAir_flow_nominal*{1, 1, 1},
      dp_nominal = {20, 20, 20},
      energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState,
      portFlowDirection_1 = Modelica.Fluid.Types.PortFlowDirection.Bidirectional,
      portFlowDirection_2 = Modelica.Fluid.Types.PortFlowDirection.Bidirectional,
      portFlowDirection_3 = Modelica.Fluid.Types.PortFlowDirection.Bidirectional,
      from_dp = true)
      annotation(Placement(transformation(origin = {-120, 96},
                                          extent = {{-12, -12}, {12, 12}})));

    Buildings.Fluid.HeatExchangers.WetCoilCounterFlow cooCoil(
      redeclare package Medium1 = MediumA,
      redeclare package Medium2 = MediumW,
      m1_flow_nominal = mAir_flow_nominal,
      m2_flow_nominal = mWat_flow_nominal,
      dp1_nominal = 150,
      dp2_nominal = 2500,
      UA_nominal = UA_nominal,
      r_nominal = r_nominal,
      nEle = nEle,
      tau1 = 3,
      tau2 = 15,
      tau_m = 8,
      energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial,
      allowFlowReversal1 = true,
      allowFlowReversal2 = true,
      from_dp1 = true,
      from_dp2 = true,
      show_T = true)
      annotation(Placement(transformation(origin = {-98, 32},
                                          extent = {{-14, -14}, {14, 14}})));

    Buildings.Fluid.Sensors.TemperatureTwoPort TsaCoo(
      redeclare package Medium = MediumA,
      m_flow_nominal = mAir_flow_nominal,
      tau = 0,
      allowFlowReversal = false)
      annotation(Placement(transformation(origin = {-62, 40},
                                          extent = {{-10, -10}, {10, 10}})));

    Buildings.Fluid.HeatExchangers.HeaterCooler_u heaCoil(
      redeclare package Medium = MediumA,
      m_flow_nominal = mAir_flow_nominal,
      dp_nominal = 75,
      Q_flow_nominal = 2000,
      energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial,
      tau = 45,
      allowFlowReversal = false,
      from_dp = false)
      annotation(Placement(transformation(origin = {-2, 40},
                                          extent = {{-12, -12}, {12, 12}})));

    Buildings.Fluid.Sensors.TemperatureTwoPort Tsa(
      redeclare package Medium = MediumA,
      m_flow_nominal = mAir_flow_nominal,
      tau = 0,
      allowFlowReversal = false)
      annotation(Placement(transformation(origin = {30, 40},
                                          extent = {{-10, -10}, {10, 10}})));

    Buildings.Fluid.Sensors.RelativeHumidityTwoPort RHsa(
      redeclare package Medium = MediumA,
      m_flow_nominal = mAir_flow_nominal,
      tau = 0,
      allowFlowReversal = false)
      annotation(Placement(transformation(origin = {58, 40},
                                          extent = {{-10, -10}, {10, 10}})));

    Buildings.Fluid.Movers.FlowControlled_m_flow fanSA(
      redeclare package Medium = MediumA,
      m_flow_nominal = mAir_flow_nominal,
      addPowerToMedium = false,
      energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState,
      allowFlowReversal = false,
      nominalValuesDefineDefaultPressureCurve = false,
      constantMassFlowRate = 0,
      dpMax = 10000,
      m_flow_start = mAir_flow_nominal*0.5,
      per(
        etaHydMet = Buildings.Fluid.Movers.BaseClasses.Types.HydraulicEfficiencyMethod.NotProvided,
        etaMotMet = Buildings.Fluid.Movers.BaseClasses.Types.MotorEfficiencyMethod.NotProvided,
        pressure(
          V_flow = {0, mAir_flow_nominal/1.2, mAir_flow_nominal/1.2*2},
          dp     = {830, 530, 0})))
      annotation(Placement(transformation(origin = {90, 40},
                                          extent = {{-12, -12}, {12, 12}})));


    Buildings.Fluid.Sensors.MassFlowRate senSA(
      redeclare package Medium = MediumA,
      allowFlowReversal = false)
      annotation(Placement(transformation(origin = {128, 40},
                                          extent = {{-10, -10}, {10, 10}})));
  

    Buildings.Fluid.MixingVolumes.MixingVolume zone(
      redeclare package Medium = MediumA,
      V = Vzone,
      nPorts = 4,
      m_flow_nominal = mAir_flow_nominal,
      energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial,
      massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial,
      p_start = 101325,
      T_start = 297.15,
      X_start = {0.008, 0.992},
      C_start = {600e-6},
      C_nominal = {1.519e-3},
      mSenFac = 10,
      allowFlowReversal = false)
      annotation(Placement(transformation(origin = {194, 54},
                                          extent = {{-14, -14}, {14, 14}})));

    Buildings.Fluid.Sensors.TemperatureTwoPort TzoneSen(
      redeclare package Medium = MediumA,
      m_flow_nominal = mAir_flow_nominal,
      tau = 0,
      allowFlowReversal = false)
      annotation(Placement(transformation(origin = {176, 8},
                                          extent = {{10, -10}, {-10, 10}})));

    Buildings.Fluid.Sensors.RelativeHumidityTwoPort RHzoneSen(
      redeclare package Medium = MediumA,
      m_flow_nominal = mAir_flow_nominal,
      tau = 0,
      allowFlowReversal = false)
      annotation(Placement(transformation(origin = {144, 8},
                                          extent = {{10, -10}, {-10, 10}})));

    Buildings.Fluid.Sensors.TraceSubstancesTwoPort CO2zoneSen(
      redeclare package Medium = MediumA,
      m_flow_nominal = mAir_flow_nominal,
      substanceName = "CO2",
      tau = 0,
      allowFlowReversal = false)
      annotation(Placement(transformation(origin = {112, 8},
                                          extent = {{10, -10}, {-10, 10}})));

    // ===== FIXED: retJ portFlowDirection_* = Bidirectional =====
    Buildings.Fluid.FixedResistances.Junction retJ(
      redeclare package Medium = MediumA,
      m_flow_nominal = mAir_flow_nominal*{1, 1, 1},
      dp_nominal = {20, 20, 20},
      energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState,
      portFlowDirection_1 = Modelica.Fluid.Types.PortFlowDirection.Bidirectional,
      portFlowDirection_2 = Modelica.Fluid.Types.PortFlowDirection.Bidirectional,
      portFlowDirection_3 = Modelica.Fluid.Types.PortFlowDirection.Bidirectional,
      linearized = true,
      from_dp = true)
      annotation(Placement(transformation(origin = {126, -42},
                                          extent = {{-12, -12}, {12, 12}})));

    Buildings.Fluid.Actuators.Dampers.Exponential damRA(
      redeclare package Medium = MediumA,
      m_flow_nominal = mAir_flow_nominal,
      dpDamper_nominal = 60,
      dpFixed_nominal = 10,
      allowFlowReversal = true,
      use_constant_density = true,
      from_dp = true)
      annotation(Placement(transformation(origin = {178, -42},
                                          extent = {{-10, -10}, {10, 10}})));

    Buildings.Fluid.Actuators.Dampers.Exponential damEA(
      redeclare package Medium = MediumA,
      m_flow_nominal = mAir_flow_nominal,
      dpDamper_nominal = 60,
      dpFixed_nominal = 10,
      allowFlowReversal = false,
      use_constant_density = true,
      from_dp = true)
      annotation(Placement(transformation(origin = {218, -92},
                                          extent = {{-10, -10}, {10, 10}})));

    Buildings.Fluid.Movers.FlowControlled_m_flow fanEA(
      redeclare package Medium = MediumA,
      m_flow_nominal = mAir_flow_nominal,
      addPowerToMedium = false,
      energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState,
      allowFlowReversal = false,
      nominalValuesDefineDefaultPressureCurve = false,
      constantMassFlowRate = 0,
      dpMax = 10000,
      m_flow_start = mAir_flow_nominal*0.3,
      per(
        etaHydMet = Buildings.Fluid.Movers.BaseClasses.Types.HydraulicEfficiencyMethod.NotProvided,
        etaMotMet = Buildings.Fluid.Movers.BaseClasses.Types.MotorEfficiencyMethod.NotProvided,
        pressure(
          V_flow = {0, mAir_flow_nominal/1.2, mAir_flow_nominal/1.2*2},
          dp     = {600, 400, 0})))
      annotation(Placement(transformation(origin = {252, -92},
                                          extent = {{-12, -12}, {12, 12}})));

    Carnot_TEva_NoWarn chiller(
      redeclare package Medium1 = MediumW,
      redeclare package Medium2 = MediumW,
      QEva_flow_nominal = QEva_flow_nominal,
      COP_nominal = COP_chiller_nominal,
      use_eta_Carnot_nominal = false,
      dp1_nominal = 3000,
      dp2_nominal = 3000,
      energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial,
      allowFlowReversal1 = false,
      allowFlowReversal2 = false,
      show_T = true,
      TEvaSet_cold = TEva_cold,
      TEvaSet_off = TEva_off)
      annotation(Placement(transformation(origin = {-180, 32},
                                          extent = {{-14, -14}, {14, 14}})));

    Buildings.Fluid.Movers.FlowControlled_dp pumpCW(
      redeclare package Medium = MediumW,
      m_flow_nominal = mWat_flow_nominal,
      nominalValuesDefineDefaultPressureCurve = true,
      energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState,
      allowFlowReversal = false,
      per(
        etaHydMet = Buildings.Fluid.Movers.BaseClasses.Types.HydraulicEfficiencyMethod.NotProvided,
        etaMotMet = Buildings.Fluid.Movers.BaseClasses.Types.MotorEfficiencyMethod.NotProvided,
        pressure(
          V_flow = {0, mWat_flow_nominal/1000},
          dp = {10000, 0})))
      annotation(Placement(transformation(origin = {-136, -62},
                                          extent = {{-12, -12}, {12, 12}})));

    Buildings.Fluid.Actuators.Valves.TwoWayLinear valCW(
      redeclare package Medium = MediumW,
      m_flow_nominal = mWat_flow_nominal,
      dpValve_nominal = 2000,
      dpFixed_nominal = 500,
      l = 0.05,
      allowFlowReversal = false,
      from_dp = false)
      annotation(Placement(transformation(origin = {-96, -62},
                                          extent = {{-10, -10}, {10, 10}})));

    Buildings.Fluid.Sources.MassFlowSource_T conSupply(
      redeclare package Medium = MediumW,
      m_flow = mWat_flow_nominal,
      T = 298.15,
      nPorts = 1)
      annotation(Placement(transformation(origin = {-260, 8},
                                          extent = {{-10, -10}, {10, 10}})));

    Buildings.Fluid.Sources.Boundary_pT conReturn(
      redeclare package Medium = MediumW,
      p = 101325,
      T = 303.15,
      nPorts = 1)
      annotation(Placement(transformation(origin = {-146, 58},
                                          extent = {{10, -10}, {-10, 10}})));

    Buildings.Fluid.Sources.Boundary_pT expVesEva(
      redeclare package Medium = MediumW,
      p = 101325,
      T = 279.15,
      nPorts = 1)
      annotation(Placement(transformation(origin = {-146, 8},
                                          extent = {{10, -10}, {-10, 10}})));

    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow intGains
      annotation(Placement(transformation(origin = {160, 84},
                                          extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Sources.RealExpression gainsSchedule(
      y = occupancy*Q_internal)
      annotation(Placement(transformation(origin = {118, 84},
                                          extent = {{-18, -12}, {18, 12}})));

    Buildings.Fluid.Sources.MassFlowSource_T CO2source(
      redeclare package Medium = MediumA,
      use_m_flow_in = true,
      use_C_in = true,
      nPorts = 1)
      annotation(Placement(transformation(origin = {242, 22},
                                          extent = {{10, -10}, {-10, 10}})));

    Modelica.Blocks.Sources.RealExpression CO2_massflow(
      y = 1e-5)
      annotation(Placement(transformation(origin = {278, 30},
                                          extent = {{10, -8}, {-10, 8}})));

    Modelica.Blocks.Sources.RealExpression CO2_concentration(
      y = max(0.0, occupancy*nOccMax*CO2_gen_per_person/1e-4))
      annotation(Placement(transformation(origin = {278, 16},
                                          extent = {{10, -8}, {-10, 8}})));

    Buildings.Fluid.Sources.MassFlowSource_T moistureSource(
      redeclare package Medium = MediumA,
      use_m_flow_in = true,
      use_X_in = true,
      nPorts = 1)
      annotation(Placement(transformation(origin = {244, 58},
                                          extent = {{10, -10}, {-10, 10}})));

    Modelica.Blocks.Sources.RealExpression moisture_gen(
      y = occupancy*nOccMax*moisture_gen_per_person)
      annotation(Placement(transformation(origin = {278, 66},
                                          extent = {{10, -8}, {-10, 8}})));

    Modelica.Blocks.Sources.RealExpression moisture_X[MediumA.nX](
      y = {1.0, 0.0})
      annotation(Placement(transformation(origin = {278, 54},
                                          extent = {{10, -6}, {-10, 6}})));

    Modelica.Blocks.Sources.RealExpression yOA_lim(
      y = max(0.2, min(1.0, uOA)))
      annotation(Placement(transformation(origin = {-263, 147},
                                          extent = {{-15, -13}, {15, 13}})));

    Modelica.Blocks.Sources.RealExpression yRA_int(
      y = max(0.05, 1 - max(0.2, min(1.0, uOA))))
      annotation(Placement(transformation(origin = {130, -20},
                                          extent = {{-12, -8}, {12, 8}})));
    Modelica.Blocks.Sources.RealExpression yEA_int(
      y = max(0.2, min(1.0, uOA)))
      annotation(Placement(transformation(origin = {158, -73},
                                          extent = {{-20, -13}, {20, 13}})));



    Modelica.Blocks.Sources.RealExpression uFan_lim(
      y = max(0.1, min(1.0, uFan))*mAir_flow_nominal)
      annotation(Placement(transformation(origin = {58, 104},
                                          extent = {{-26, -12}, {26, 12}})));

    Modelica.Blocks.Sources.RealExpression uFanEA_lim(
      y = max(0.05, min(1.0, uFanEA))*mAir_flow_nominal)
      annotation(Placement(transformation(origin = {226, -38},
                                          extent = {{-16, -12}, {16, 12}})));



    Modelica.Blocks.Sources.Constant dpPump_set(k = 8000)
      annotation(Placement(transformation(origin = {-176, -22},
                                          extent = {{-10, -10}, {10, 10}})));

    Buildings.Controls.Continuous.LimPID pidCW(
      controllerType = Modelica.Blocks.Types.SimpleController.PI,
      k = 0.3,
      Ti = 300,
      yMax = 0.95,
      yMin = 0.05,
      y_start = 0.5,
      xi_start = 0.5,
      reverseActing = true,
      strict = false,
      initType = Modelica.Blocks.Types.Init.InitialOutput)
      annotation(Placement(transformation(origin = {-62, 72},
                                          extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Sources.Constant Tsa_cold_set(k = 285.15)
      annotation(Placement(transformation(origin = {-92, 72},
                                          extent = {{-10, -10}, {10, 10}})));

    Modelica.Blocks.Math.Product heaterEnable
      annotation(Placement(transformation(origin = {-34, 120},
                                          extent = {{-6, -6}, {6, 6}})));

    Modelica.Blocks.Sources.RealExpression flowFraction(
      y = min(1.0, max(0.0, senSA.m_flow/(mAir_flow_nominal*0.15))))
      annotation(Placement(transformation(origin = {-85, 116},
                                          extent = {{-15, -14}, {15, 14}})));

    Modelica.Blocks.Sources.RealExpression uHeater_lim(
      y = max(0.0, min(1.0, uHeater)))
      annotation(Placement(transformation(origin = {-85, 136},
                                          extent = {{-15, -14}, {15, 14}})));

    Modelica.Blocks.Sources.Constant TBase_eva(k = T_coil_dehum)
      annotation(Placement(transformation(origin = {-269.333, 44.8},
                                          extent = {{-21.3333, -12.8}, {21.3333, 12.8}})));

    Modelica.Blocks.Sources.RealExpression uChiller_lim(
      y = max(0.0, min(1.0, uChiller)))
      annotation(Placement(transformation(origin = {-263, -18},
                                          extent = {{-15, -12}, {15, 12}})));

    Real rhoSA(nominal = 1.2);

  equation
    // *** CONTROL SIGNALS - DARK BLUE {0,0,127} ***
    connect(TDryBul_in, weaDat.TDryBul_in) annotation(
      Line(points = {{-366, 122}, {-335, 122}, {-335, 121}},
           color = {0, 0, 127}));
    connect(relHum_in, weaDat.relHum_in) annotation(
      Line(points = {{-366, 108}, {-365.5, 108}, {-365.5, 110}, {-335, 110}},
           color = {0, 0, 127}));
    connect(pAtm_in, weaDat.pAtm_in) annotation(
      Line(points = {{-366, 136}, {-365, 136}, {-365, 134}, {-335, 134}},
           color = {0, 0, 127}));
    connect(HGloHor_in, weaDat.HGloHor_in) annotation(
      Line(points = {{-366, 58}, {-365, 58}, {-365, 60}, {-335, 60}},
           color = {0, 0, 127}));
    connect(HDifHor_in, weaDat.HDifHor_in) annotation(
      Line(points = {{-366, 72}, {-365, 72}, {-365, 69}, {-335, 69}},
           color = {0, 0, 127}));
    connect(winSpe_in, weaDat.winSpe_in) annotation(
      Line(points = {{-366, 86}, {-365, 86}, {-365, 85}, {-335, 85}},
           color = {0, 0, 127}));

    // *** WEATHER BUS - YELLOW/ORANGE {255,204,51} ***
    connect(weaDat.weaBus, ambAir.weaBus) annotation(
      Line(points = {{-276, 96}, {-260, 96}},
           color = {255, 204, 51},
           thickness = 0.5));

    // ============ CONTROL CONNECTIONS ============
    connect(yOA_lim.y, filtOA.u) annotation(
      Line(points = {{-250, 133}, {-234, 133}},
           color = {0, 0, 127}));
    connect(filtOA.y, limOA.u) annotation(
      Line(points = {{-222, 133}, {-214, 133}},
           color = {0, 0, 127}));
    connect(limOA.y, damOA.y) annotation(
      Line(points = {{-202, 133}, {-210, 133}, {-210, 92}, {-192, 92}},
           color = {0, 0, 127}));

    connect(yRA_int.y, filtRA.u) annotation(
      Line(points = {{143.2, -20}, {155.2, -20}},
           color = {0, 0, 127}));
    connect(filtRA.y, damRA.y) annotation(
      Line(points = {{168.6, -20}, {177.6, -20}, {177.6, -30}},
           color = {0, 0, 127}));

    connect(yEA_int.y, filtEA.u) annotation(
      Line(points = {{176, -73}, {192, -73}},
           color = {0, 0, 127}));
    connect(filtEA.y, limEA.u) annotation(
      Line(points = {{204, -73}, {212, -73}},
           color = {0, 0, 127}));
    connect(limEA.y, damEA.y) annotation(
      Line(points = {{224, -73}, {214, -73}},
           color = {0, 0, 127}));

    connect(Tsa_cold_set.y, pidCW.u_s) annotation(
      Line(points = {{-81, 72}, {-74, 72}},
           color = {0, 0, 127}));
    connect(TsaCoo.T, pidCW.u_m) annotation(
      Line(points = {{-62, 51}, {-62, 60}},
           color = {0, 0, 127}));
    connect(pidCW.y, filtValve.u) annotation(
      Line(points = {{-51, 72}, {-48, 72}, {-48, 62}, {-45, 62}},
           color = {0, 0, 127}));
    connect(filtValve.y, valCW.y) annotation(
      Line(points = {{-31.4, 62}, {-28.4, 62}, {-28.4, -22}, {-95.8, -22}, {-95.8, -50}, {-96.8, -50}},
           color = {0, 0, 127}));

    connect(uHeater_lim.y, filtHeater.u) annotation(
      Line(points = {{-68.5, 136}, {-63, 136}},
           color = {0, 0, 127}));
    connect(filtHeater.y, heaterEnable.u1) annotation(
      Line(points = {{-49.4, 136}, {-47, 136}, {-47, 124}, {-41.4, 124}},
           color = {0, 0, 127}));
    connect(flowFraction.y, heaterEnable.u2) annotation(
      Line(points = {{-68.5, 116}, {-41, 116}},
           color = {0, 0, 127}));
    connect(heaterEnable.y, heaCoil.u) annotation(
      Line(points = {{-27.4, 120}, {-27.4, 119}, {-20.3, 119}, {-20.3, 47}, {-16.4, 47}},
           color = {0, 0, 127}));

    connect(uFan_lim.y, fanSA.m_flow_in) annotation(
      Line(points = {{86.6, 104}, {90.2, 104}, {90.2, 54}, {89.2, 54}},
           color = {0, 0, 127}));

    connect(uFanEA_lim.y, fanEA.m_flow_in) annotation(
      Line(points = {{243.6, -38}, {252, -38}, {252, -78}},
           color = {0, 0, 127}));



    connect(dpPump_set.y, pumpCW.dp_in) annotation(
      Line(points = {{-165, -22}, {-136, -22}, {-136, -48}},
           color = {0, 0, 127}));

    connect(TBase_eva.y, chiller.TSet) annotation(
      Line(points = {{-245.866, 44.8}, {-196.866, 44.8}},
           color = {0, 0, 127}));
    connect(uChiller_lim.y, chiller.u) annotation(
      Line(points = {{-246.5, -18}, {-206, -18}, {-206, 40}, {-197, 40}},
           color = {0, 0, 127}));

    connect(gainsSchedule.y, intGains.Q_flow) annotation(
      Line(points = {{137.8, 84}, {149.8, 84}},
           color = {0, 0, 127}));

    connect(CO2_massflow.y, CO2source.m_flow_in) annotation(
      Line(points = {{267, 30}, {254, 30}},
           color = {0, 0, 127}));
    connect(CO2_concentration.y, CO2source.C_in[1]) annotation(
      Line(points = {{267, 16}, {254, 16}},
           color = {0, 0, 127}));

    connect(moisture_gen.y, moistureSource.m_flow_in) annotation(
      Line(points = {{267, 66}, {256, 66}},
           color = {0, 0, 127}));
    connect(moisture_X.y, moistureSource.X_in) annotation(
      Line(points = {{267, 54}, {256, 54}},
           color = {0, 0, 127}));

    // ============ FLUID CONNECTIONS ============
    connect(ambAir.ports[1], damOA.port_a) annotation(
      Line(points = {{-236, 96}, {-226, 96}},
           color = {221, 203, 0},
           thickness = 1));
    connect(damOA.port_b, senOA.port_a) annotation(
      Line(points = {{-206, 96}, {-196, 96}},
           color = {221, 203, 0},
           thickness = 1));
    connect(senOA.port_b, preFil.port_a) annotation(
      Line(points = {{-176, 96}, {-166, 96}},
           color = {221, 203, 0},
           thickness = 1));
    connect(preFil.port_b, mixT.port_1) annotation(
      Line(points = {{-146, 96}, {-132, 96}},
           color = {221, 203, 0},
           thickness = 1));
    connect(mixT.port_3, cooCoil.port_b1) annotation(
      Line(points = {{-120, 84}, {-120, 40}, {-112, 40}},
           color = {85, 255, 0},
           thickness = 1));
    connect(cooCoil.port_a1, TsaCoo.port_a) annotation(
      Line(points = {{-84, 40.4}, {-72, 40.4}},
           color = {85, 255, 0},
           thickness = 1));
    connect(TsaCoo.port_b, heaCoil.port_a) annotation(
      Line(points = {{-52, 40}, {-14, 40}},
           color = {85, 255, 0},
           thickness = 1));
    connect(heaCoil.port_b, Tsa.port_a) annotation(
      Line(points = {{10, 40}, {20, 40}},
           color = {85, 255, 0},
           thickness = 1));
    connect(Tsa.port_b, RHsa.port_a) annotation(
      Line(points = {{40, 40}, {48, 40}},
           color = {0, 255, 0},
           thickness = 1));
    connect(RHsa.port_b, fanSA.port_a) annotation(
      Line(points = {{68, 40}, {78, 40}},
           color = {85, 255, 0},
           thickness = 1));
    connect(fanSA.port_b, senSA.port_a) annotation(
        Line(points = {{102, 40}, {118, 40}},
             color = {85, 255, 0},
             thickness = 1));
    connect(senSA.port_b, zone.ports[1]) annotation(
        Line(points = {{138, 40}, {194, 40}},
             color = {85, 255, 0},
             thickness = 1));
  
    connect(zone.ports[2], TzoneSen.port_a) annotation(
      Line(points = {{194, 40}, {194, 8}, {186, 8}},
           color = {249, 171, 255},
           thickness = 1));
    connect(TzoneSen.port_b, RHzoneSen.port_a) annotation(
      Line(points = {{166, 8}, {154, 8}},
           color = {249, 171, 255},
           thickness = 1));
    connect(RHzoneSen.port_b, CO2zoneSen.port_a) annotation(
      Line(points = {{134, 8}, {122, 8}},
           color = {249, 171, 255},
           thickness = 1));
    connect(CO2zoneSen.port_b, retJ.port_1) annotation(
      Line(points = {{102, 8}, {88, 8}, {88, -42}, {114, -42}},
           color = {249, 171, 255},
           thickness = 1));
    connect(retJ.port_2, damRA.port_b) annotation(
      Line(points = {{138, -42}, {168, -42}},
           color = {255, 157, 0},
           thickness = 1));
    connect(damRA.port_a, mixT.port_2) annotation(
      Line(points = {{188, -42}, {218, -42}, {218, 96}, {-108, 96}},
           color = {255, 157, 0},
           thickness = 1));
    connect(retJ.port_3, damEA.port_a) annotation(
      Line(points = {{126, -54}, {126, -92}, {208, -92}},
           color = {255, 85, 0},
           thickness = 1));
    connect(damEA.port_b, fanEA.port_a) annotation(
      Line(points = {{228, -92}, {240, -92}},
           color = {255, 85, 0},
           thickness = 1));
    connect(fanEA.port_b, ambAir.ports[2]) annotation(
      Line(points = {{264, -92}, {270, -92}, {270, -112}, {-236, -112}, {-236, 96}},
           color = {255, 85, 0},
           thickness = 1));
    connect(CO2source.ports[1], zone.ports[3]) annotation(
      Line(points = {{232, 22}, {232, 21.5}, {226, 21.5}, {226, 39}, {194, 39}, {194, 40}},
           color = {0, 127, 255},
           thickness = 1));
    connect(moistureSource.ports[1], zone.ports[4]) annotation(
      Line(points = {{234, 58}, {234, 58.5}, {226, 58.5}, {226, 39}, {194, 39}, {194, 40}},
           color = {0, 127, 255},
           thickness = 1));

    connect(conSupply.ports[1], chiller.port_a1) annotation(
      Line(points = {{-250, 8}, {-221, 8}, {-221, 40}, {-194, 40}},
           color = {0, 127, 255},
           thickness = 1));
    connect(chiller.port_b1, conReturn.ports[1]) annotation(
      Line(points = {{-166, 40.4}, {-161, 40.4}, {-161, 58.4}, {-156, 58.4}},
           color = {0, 127, 255},
           thickness = 1));
    connect(expVesEva.ports[1], chiller.port_a2) annotation(
      Line(points = {{-156, 8}, {-166, 8}, {-166, 24}},
           color = {170, 255, 255},
           thickness = 1));
    connect(chiller.port_b2, pumpCW.port_a) annotation(
      Line(points = {{-194, 23.6}, {-199, 23.6}, {-199, -62.4}, {-148, -62.4}},
           color = {85, 255, 255},
           thickness = 1));
    connect(pumpCW.port_b, valCW.port_a) annotation(
      Line(points = {{-124, -62}, {-106, -62}},
           color = {85, 255, 255},
           thickness = 1));
    connect(valCW.port_b, cooCoil.port_a2) annotation(
      Line(points = {{-86, -62}, {-72, -62}, {-72, 24}, {-84, 24}},
           color = {85, 255, 255},
           thickness = 1));
    connect(cooCoil.port_b2, chiller.port_a2) annotation(
      Line(points = {{-112, 23.6}, {-166, 23.6}},
           color = {85, 255, 255},
           thickness = 1));
    connect(intGains.port, zone.heatPort) annotation(
      Line(points = {{170, 84}, {180, 84}, {180, 54}},
           color = {191, 0, 0}));

    // ============ OUTPUT CALCULATIONS ============
    rhoSA = max(
      0.5,
      min(
        2.0,
        Buildings.Media.Air.density(
          Buildings.Media.Air.setState_pTX(
            101325,
            max(273.15, min(323.15, Tsa.T)),
            {0.01, 0.99}))));

    T_SA = Tsa.T;
    RH_SA = max(0.0, min(1.0, RHsa.phi));
    Vdot_SA = max(0.0, senSA.m_flow)/rhoSA;

    T_zone = TzoneSen.T;
    RH_zone = max(0.0, min(1.0, RHzoneSen.phi));
    CO2_zone_ppm = max(0.0, CO2zoneSen.C)*1e6;

    P_fan   = Pfan_nominal  * max(0.0, min(1.0, fanSA.m_flow/mAir_flow_nominal))^3;
    P_fanEA = PfanEA_nominal* max(0.0, min(1.0, fanEA.m_flow/mAir_flow_nominal))^3;
    Q_chiller = max(0.0, -chiller.QEva_flow);
    P_chiller = max(0.0, chiller.P);
    P_pump = PpumpCW_nominal*max(0.0, min(1.0, pumpCW.m_flow/mWat_flow_nominal))^3;

    Q_heater = heaCoil.Q_flow;
    T_SA_afterCooling = TsaCoo.T;

    P_total = P_fan + P_fanEA + P_chiller + P_pump;

    annotation(
      experiment(
        StartTime = 0,
        StopTime = 31536000,
        Interval = 900,
        Tolerance = 1e-06),
      Diagram(coordinateSystem(extent = {{-380, 180}, {340, -120}})),
      Documentation(info = "<html>
<p><b>HVAC MODEL - RL-READY with External Weather Data Inputs</b></p>
<p>Model accepts weather data from external inputs for Python-driven simulation.</p>
</html>"));
  end AHU_FMU_Core;

  annotation(
    uses(
      Modelica(version = "4.0.0"),
      Buildings(version = "12.1.0")),
    Documentation(info = "<html>
<p><b>HVAC_FMU Package - RL-READY with External Weather Data</b></p>
<p>Model receives weather data through RealInput connectors for external control, outdoor CO2 fixed at 400 ppm.</p>
</html>"));
end HVAC_FMU;
