Parameters described on this page:
[SRC_DELEPTONIZATION](#SRC_DELEPTONIZATION), &nbsp;
[SRC_TURBULENCE](#SRC_TURBULENCE), &nbsp;
[SRC_USER](#SRC_USER), &nbsp;
[TURB_VEL](#TURB_VEL), &nbsp;
[TURB_AMPL_FACTOR](#TURB_AMPL_FACTOR), &nbsp;
[TURB_KDRIV](#TURB_KDRIV), &nbsp;
[TURB_KMIN](#TURB_KMIN), &nbsp;
[TURB_KMAX](#TURB_KMAX), &nbsp;
[TURB_ZETA](#TURB_ZETA), &nbsp;
[TURB_SPEC_FORM](#TURB_SPEC_FORM), &nbsp;
[TURB_POW](#TURB_POW), &nbsp;
[TURB_RSEED_INIT](#TURB_RSEED_INIT), &nbsp;
[TURB_UPDATE_STEP](#TURB_UPDATE_STEP), &nbsp;
[TURB_TABLE_SIZE](#TURB_TABLE_SIZE) &nbsp;
[TURB_RESET](#TURB_RESET), &nbsp;

Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="SRC_DELEPTONIZATION"></a>
* #### `SRC_DELEPTONIZATION` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Deleptonization (for simulations of stellar core collapse).
    * **Restriction:**
Only applicable when enabling the compilation option
[[--model | [Installation]-Option-List#--model]]=HYDRO.

<a name="SRC_TURBULENCE"></a>
* #### `SRC_TURBULENCE` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Turbulence acceleration source term.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--model | [Installation]-Option-List#--model]]=HYDRO.

<a name="SRC_USER"></a>
* #### `SRC_USER` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
User-defined source terms.
See [[Add User-defined Source Terms | Source-Terms#add-user-defined-source-terms]] for details.
    * **Restriction:**

<a name="TURB_VEL"></a>
* #### `TURB_VEL` &ensp; (>0.0) &ensp; [0.2]
    * **Description:**
Target turbulence velocity dispersion in code unit.
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="TURB_AMPL_FACTOR"></a>
* #### `TURB_AMPL_FACTOR` &ensp; (>0.0) &ensp; [1.0]
    * **Description:**
Amplitude multiplier for turbulence acceleration field.
Adjust this factor if velocity dispersion measured is different from target velocity dispersion.
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="TURB_KDRIV"></a>
* #### `TURB_KDRIV` &ensp; (>0.0) &ensp; [2.0]
    * **Description:**
Turbulence driving wave number in units of 2π/BOX_SIZE, determine the correlation time τ=BOX_SIZE/(TURB_KDRIV*TURB_VEL)
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="TURB_KMIN"></a>
* #### `TURB_KMIN` &ensp; (&#8805;1.0) &ensp; [1.0]
    * **Description:**
Turbulence minimum wave number in units of 2π/BOX_SIZE.
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="TURB_KMAX"></a>
* #### `TURB_KMAX` &ensp; (&#8805;1.0) &ensp; [3.0]
    * **Description:**
Turbulence maximum wave number in units of 2*$\pi$/BOX_SIZE.
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="TURB_ZETA"></a>
* #### `TURB_ZETA` &ensp; (1.0&#8805;input&#8805;0.0) &ensp; [1.0]
    * **Description:**
Turbulence solenoidal weighting (0.0=fully compressive, 1.0=fully solenoidal).
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="TURB_SPEC_FORM"></a>
* #### `TURB_SPEC_FORM` &ensp; (0=constant, 1=parabolic, 2=power law) &ensp; [1]
    * **Description:**
Spectral form of the turbulence driving amplitude.
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="TURB_POW"></a>
* #### `TURB_POW` &ensp; (none) &ensp; [-5/3]
    * **Description:**
Power law for energy power spectrum.
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1, [TURB_SPEC_FORM](#SRC_TURB_SPEC_FORM)=2.

<a name="TURB_RSEED_INIT"></a>
* #### `TURB_RSEED_INIT` &ensp; (>0) &ensp; [123]
    * **Description:**
Initial random seed for turbulence, useless when restart without reset.
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="TURB_UPDATE_STEP"></a>
* #### `TURB_UPDATE_STEP` &ensp; (&#8805;1) &ensp; [10]
    * **Description:**
Update turbulence pattern every dt=(correlation time/TURB_UPDATE_STEP).
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="TURB_TABLE_SIZE"></a>
* #### `TURB_TABLE_SIZE` &ensp; (&#8805;0) &ensp; [128]
    * **Description:**
Resolution of turbulence table.
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

<a name="TURB_RESET"></a>
* #### `TURB_RESET` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Reset turbulence field when restart.
    * **Restriction:**
Only applicable when [[--model | [Installation]-Option-List#--model]]=HYDRO and
[SRC_TURBULENCE](#SRC_TURBULENCE)=1.

## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Main page of Source Terms | Source-Terms]]
