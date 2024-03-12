'''The Purpose of this file is simply to take data from the SQUID-data-sample-bank
google doc and write it to a pandas dataframe'''


import gspread
import pandas as pd
from oauth2client.service_account import ServiceAccountCredentials
#Code to authorize access to google spreadsheet
scope = ['https://spreadsheets.google.com/feeds','https://www.googleapis.com/auth/drive']
creds = ServiceAccountCredentials.from_json_keyfile_name('testsheetkey.json', scope)
client = gspread.authorize(creds)

#Code to retrieve data from google sheet via gspread library
sheet = client.open('STFMR-data-sample-bank').sheet1
values = sheet.get_all_values()

#Code to transform data from gspread into pandas table
raw_table = pd.DataFrame(values[1:],columns=values[0])
#Code to Set the sample directory as the index for the pandas dataframe
data = raw_table.set_index('Directory')

print(data)
#everything below was for testing purposes - but I will leave it for the syntax

#print(data)
#fig, ax = plt.subplots()
#ax.plot([np.float(x) for x in data['Mass (g)']])
#plt.show()
#print(data.index)
'''filename = 'LSMO_071321a_IP.dat'



temp = data.loc[filename]
print(temp)'''
testsheet = sheet.get_all_records() #This just shows what gspread outputs the data as, whereas
#it is easier to use pandas
#print(testsheet)