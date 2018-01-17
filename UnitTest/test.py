from FRUIT import *

#test_modules = ["Library/allp_mod_test.f90 Library/all_mod_test.f90"]
test_modules = ["Library/allp_mod_test.f90", "Library/all_mod_test.f90", 				"Process/Calc_Sgs_Coefficients_Dynamic_Smagorinsky_test.f90"]
driver = "test_driver.f90"
#build_command="make test_driver"

# run_command="./allp_mod"

suite = test_suite(test_modules)

# build test
#suite.build(driver, build_command)

suite.build_run(driver)

# run test
#suite.run(driver, run_command)

suite.summary()

