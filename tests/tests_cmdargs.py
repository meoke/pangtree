import unittest
import argparse
from ddt import ddt, data, unpack

from context import cmdargs


@ddt
class CmdargsTest(unittest.TestCase):

    @data("", "foo")
    def test_file_arg_incorrect_path(self, value):
        with self.assertRaises(argparse.ArgumentTypeError) as err:
            actual_output = cmdargs._file_arg(value)
            self.assertEqual(str(err.exception), f"File {value} does not exist or is not a file.")

    @data("1", "0", "0.01")
    def test_float_0_1_type_correct_value(self, value):
        actual_output = cmdargs._float_0_1(value)
        self.assertEqual(actual_output, float(value))

    @data("-1", "-1.5", "1.01", "10")
    def test_float_0_1_type_value_out_of_range(self, value):
        with self.assertRaises(argparse.ArgumentTypeError) as err:
            actual_output = cmdargs._float_0_1(value)
            self.assertEqual(str(err), 'This argument must be in range [0,1].')

    @data("A", "1.5.6", "", "10")
    def test_float_0_1_type_not_float(self, value):
        with self.assertRaises(argparse.ArgumentTypeError) as err:
            actual_output = cmdargs._float_0_1(value)
            self.assertEqual(str(err), f"{value} was passed, a float excpected.")


if __name__ == '__main__':
    unittest.main()
