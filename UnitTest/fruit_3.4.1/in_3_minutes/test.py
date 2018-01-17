from FRUIT import *

test_modules = ['calculator_test.f90']
driver = "calculator_test_driver.f90"

suite = test_suite(test_modules)
suite.build_run(driver)
suite.summary()
