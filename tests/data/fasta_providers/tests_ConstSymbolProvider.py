import unittest

from ...context import ConstSymbolProvider, MissingSymbol, InputError
from ...context import pSeq, pNode

class ConstSymbolProviderTests(unittest.TestCase):

    def test_1_no_symbol_provided(self):
        missing_symbol = MissingSymbol()
        const_symbol_provider = ConstSymbolProvider(missing_symbol)

        expected_symbol = pNode.Base('?')
        actual_symbol = const_symbol_provider.get_base(pSeq.SequenceID('s'), 0)
        self.assertEqual(expected_symbol, actual_symbol)

    def test_2_symbol_provided(self):
        const_symbol_provider = ConstSymbolProvider(MissingSymbol('*'))

        expected_symbol =pNode.Base('*')
        actual_symbol = const_symbol_provider.get_base(pSeq.SequenceID('s'), 0)
        self.assertEqual(expected_symbol, actual_symbol)

    def test_3_incorrect_missing_symbol(self):
        with self.assertRaises(InputError) as e:
            _ = MissingSymbol('**')

        expected_message = 'Missing symbol must be a single character.'
        actual_message = str(e.exception)
        self.assertEqual(expected_message, actual_message)


if __name__ == '__main__':
    unittest.main()
