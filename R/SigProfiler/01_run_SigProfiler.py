def main_function():
	# import sys
	# print(sys.path)

	import argparse
	from SigProfilerExtractor import sigpro as sig

	## Args
	parser = argparse.ArgumentParser()
	parser.add_argument("--input_data", dest="input_data", type=str, help="Path to mutation context matrix txt")
	parser.add_argument("--context_type", dest="context_type", type=str, default="96", help="Context type (can be '96','DINUC','ID')")
	parser.add_argument("--output", dest="output", type=str, help="Output directory")
	parser.add_argument("--minimum_signatures", dest="minimum_signatures", type=int, default=1, help="Minimum number of signatures")
	parser.add_argument("--maximum_signatures", dest="maximum_signatures", type=int, default=21, help="Maximum number of signatures")
	args = parser.parse_args()


	## Main
	sig.sigProfilerExtractor(
		output=args.output, 
		input_data=args.input_data,
		input_type="matrix",
		context_type=args.context_type,
		reference_genome="GRCh37", 
		minimum_signatures=args.minimum_signatures, 
		maximum_signatures=args.maximum_signatures, 
		nmf_replicates=100, 
		cpu=-1
	)

if __name__=="__main__":
	main_function()
