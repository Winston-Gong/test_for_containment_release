# Contact network based simulation of COVID-19 outbreaks

Examples that illustrate and code that is generated from the study "**Modelling the impact of rapid tests, tracing and distancing in lower-income countries suggest optimal policies varies with rural-urban settings**": https://www.medrxiv.org/content/10.1101/2021.03.17.21253853v1 

## Quick example

Example code in local_example.R provide quick local simulation for demonstration. Run the script to generate a synthetic population and outbreaks on it. 

Example simulation have following specification:

1. Population size = 200

2. Default community setting = rural

3. Number of outbreaks = 20

Running the example file local_example.R  will infer one small community and simulate outbreaks on it. Two graphs will be saved: one represent an example of data contact structure; the other represent the average daily contact number in this community. 

## Code structure

The simulation contain two steps: linux_generate_network.R will generate several synthetic populations and infer the contact network for each population; linux_simulate.R will simulate outbreaks and containment over the synthetic populations. Supplementary Table 4 contains descriptions of all simulation settings considered in the study, which are listed in linux_simulate.R.

Multisetting_20200519_base.R simulate 200 bootstrap trajectories for four settings: baseline (No NPI imposed), type 1 combined strategy (symptom-based diagnostic + PCR), type 2 combined strategy (symptom-based diagnostic + PCR + antibody RDT) and type 3 combined strategy (antigen RDT  + PCR ). Functions that simulated transmission of virus, daily contact, testing and quarantine and saved in Multisetting_20200519_functions.R.  

For researcher who aim to adapt the code for their own community simulation, Table 1 and Table 2 in the paper contains reference of our parameter settings, which are coded in linux_generate_network.R and linux_simulate.R in using readable varialbe names. 

* Note: We recommend using the demographic paramenter settings for the three community settings. The demographic parameters (contact number, household sizes, age mixing) are highly correlated and changing those parameters might prevent the ERGM to converge. If researchers prefer to simulate using their own community setting, we strongly recommend starting with our parameter setting and change only one parameter at a time. To illustrate an example of community that is not realistic, one could imaging if we set a small household sizes and prevent out-of-household contact simultaneously. In this case it would not be possible to have high contact rate among kids since there are not enough kids within each household for contact. When simulating reduced out-of-household contact, we recomment using the social_distancing_flg parameter as this will scale the contact age mixing as well. 

## Full scale simulation on Linux-based cluster

To replicate the simulations in the study require relatively abundant compuation resources. Reaserchers should create two directories in the repository: Networks/ and Dynamics/ ; the simulation should be performed using current directory.

Code in linux_generate_network.R should be run before linux_simulate.R to create the synthetic population and infer contact structure, which will be saved in Networks/; once this step is finished, linux_simulate.R will simulate outbreaks under different containment strategy. 

Researchers aim to adapt the code should be able to submit the code to their own computation facillates. Five handles are provided for external output which should be obvious for those aiming to run their own analysis:

* Community setting: 1 for rural; 2 for non-slum urban; 3 for slum

* v_id: Choosing one varialbe to alter. Refer to Supplementary Table 4 for all the availabe choices. 

* s_id: Choosing a value for the v_id. Refer to Supplementary Table 4 for all the availabe choices.

* social_distancing_flg: choosing how much physical distancing to impose (available choices: 1,2,3,4,5); 1 for no physical distancint; 2-5 will reduce out-of-household contact number with 20% increment each.

* country ID: choosing which contry to use as reference for demographic information: 1 for Uganda, 2 for South africa, 3 for Kenya and 4 for Nigeria.
