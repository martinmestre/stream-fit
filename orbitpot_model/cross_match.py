# Perform cross match between Malhan+19 table and Gaia DR2 data
# downloaded source by source from Gaia archive.

def cross_match(table_file):
	from astropy.io import ascii
	data = ascii.read(table_file)  
	print(type(data))

	arr = data['Pmember']
	boolarr=  (arr=='Y')

	# Now we work with Gaia archive
	from astroquery.gaia import Gaia
	gaiadr2_table = Gaia.load_table('gaiadr2.gaia_source')
	for column in (gaiadr2_table.columns):
	  print(column.name,column.unit)

	# We create new empty columns to the previous table
	data.add_column(0.0, name='pmra')  # Add column of all 1.0 to end of table
	data.add_column(0.0, name='pmra_error')
	data.add_column(0.0, name='pmdec')
	data.add_column(0.0, name='pmdec_error')

	radius = u.Quantity(0.0001, u.deg)
	for i in range(0,len(data['RAJ2000'])):
		coords = SkyCoord(ra=data['RAJ2000'][i], dec=data['DEJ2000'][i], unit=(u.degree, u.degree), frame='icrs')
		j = Gaia.cone_search_async(coords, radius)
		g = j.get_results()
		data['pmra'][i]=g['pmra']
		data['pmra_error'][i]=g['pmra_error']
		data['pmdec'][i]=g['pmdec']
		data['pmdec_error'][i]=g['pmdec_error']

	gdr2_ra=(data[boolarr])['RAJ2000']
	gdr2_dec=(data[boolarr])['DEJ2000']
	gdr2_pmra=(data[boolarr])['pmra']
	gdr2_pmdec=(data[boolarr])['pmdec']

	# Save cross-matched data
	ascii.write(data, 'cross_matched.dat', overwrite=True)
	return data
