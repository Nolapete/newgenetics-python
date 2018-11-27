#just a holder for the traits and genotypes of a parent
#the constructor can be reformated, as long as self.types holds a dictionary of trait names to genotype tuples
#example entry: {"Picasso": ('P','+')}
class parent :
	#example with current input: Picasso would be ('Picasso', ('P', '+'))
	def __init__ (self, types) :
		self.types = {}
		for t in types :
			self.types[t[0]] = t[1]

	#return our type dictionary
	def get_types(self) :
		print("Types_")
		print(self.types)
		return self.types


#given p1 and p2, determine which one should be in front, considering wildtypes and recessive vs dominant
#any genotype that starts with a lowercase letter is considered recssive for this
def format_allele(p1, p2) :
	if(p1 == '+') :
		if(ord(p2[0]) >= ord('a') and ord(p2[0]) <= ord('z')) :
			return p1 + '/' + p2
		return p2 + '/' + p1
	elif (p2 == '+') :
		if(ord(p1[0]) >= ord('a') and ord(p1[0]) <= ord('z')) :
			return p2 + '/' + p1
	return p1 + '/' + p2

#uses the properties of an expanded dihybrid cross table to determine the resultant type
#ind is a tuple of (x, y) for horizontal and vertical index
#length is h where the table is h x h
#h_axis is the genotypes along the horizontal axis
#v_axis is the genotypes along the vertical axis
#h_axis and v_axis should be sorted in the same order
def cross_at_index(ind, length, h_axis, v_axis) :
	#start_x and start_y are used to determine which genotype to take at a given index
	start_x = length // 2
	start_y = length // 2

	#resulting genotype string
	res = ''

	x = ind[0]
	y = ind[1]

	#the x and y components of the genotype as calculated by the whileloop
	xpart = None
	ypart = None
	#our current index in the axes
	ax_ind = 0
	#while we still have elements left, determine which component to take based on the index
	while(start_x > 0) :
		if(x >= start_x) :
			xpart = h_axis[ax_ind][1]
		else :
			xpart = h_axis[ax_ind][0]
		if(y >= start_y) :
			ypart = v_axis[ax_ind][1]
		else :
			ypart = v_axis[ax_ind][0]

		#cut out full wildtypes
		if(xpart != '+' or ypart != '+') :
			res += ' ' + format_allele(xpart, ypart)

		#setup for the next loop
		x %= start_x
		y %= start_y
		start_x //= 2
		start_y //= 2
		ax_ind += 1

	#clear extra whitespace and return the genotype
	return res.strip()


#f1 and f2 are parent objects
def cross_parent(f1, f2) :
	f1_types = f1.get_types()
	f1_type_names = f1_types.keys()
	f2_types = f2.get_types()
	f2_type_names = f2_types.keys()

	#find any missing keys from either parent
	f2_missing = list(filter(lambda a: not a in f2_type_names, f1_type_names))
	f1_missing = list(filter(lambda a: not a in f1_type_names, f2_type_names))

	#add in the missing keys
	for n in f1_missing :
		f1_types[n] = ('+','+')

	for n in f2_missing :
		f2_types[n] = ('+','+')

	f1_final = []
	f2_final = []
	#based on the ordering of one of the parent dictionaries, create a horizontal and vertical axis list
	#the lists need to be in the same order, so using one dictionary's ordering ensures this
	names = list(f1_types.keys())
	for n in names :
		f1_final.append(f1_types[n])
		f2_final.append(f2_types[n])

	#our table grows exponentially with the number of traits to consider
	table_length = 2 ** len(f1_final)

	#for tracking the percentages for each genotype
	total_entries = float(table_length * table_length)
	percent_counter = {}

	#loop through the entire cross table, adding up the percentages
	for x in range(0, table_length) :
		for y in range(0, table_length) :
			res = cross_at_index((x, y), table_length, f1_final, f2_final)
			if (res in percent_counter) :
				percent_counter[res] += 100. / total_entries
			else :
				percent_counter[res] = 100. / total_entries

	return percent_counter



if __name__ == '__main__' :
	f1 = parent([('Sm', ('Sm', '+')), ('Z', ('Z', '+')), ('St', ('St','+'))])
	f2 = parent([('S', ('S', '+')), ('pb', ('pb', 'pb')), ('p', ('p', 'p'))])
	results = cross_parent(f1, f2)

	for key in results.keys() :
		print (str(key) + ' : ' + str(results[key]))
