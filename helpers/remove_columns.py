import pandas as pd

# Load the metadata and filtered data
metadata = pd.read_csv("static/examples/metadata_myeloma.csv")
filtered_data = pd.read_csv("filtered_data.csv")

# Extract GSM identifiers from the filtered data
# Assuming the first row in filtered_data.csv are the GSMs (columns after 'ID_REF')
filtered_gsms = filtered_data.columns[1:]

# Filter the metadata to keep only the rows where GSM is in the filtered GSMs
filtered_metadata = metadata.loc[:, ["Sample_geo_accession"] + list(filtered_gsms)]

# Save the filtered metadata to a new CSV file
filtered_metadata.to_csv("filtered_metadata.csv", index=False)

print(f"Filtered metadata saved to 'filtered_metadata.csv' with {len(filtered_metadata.columns) - 1} GSMs.")
