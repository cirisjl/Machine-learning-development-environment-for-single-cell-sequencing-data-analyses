// init-mongo.js

// Connect to MongoDB
conn = new Mongo();
db = conn.getDB("oscb");

// Check if documents already exist in the collection
var countDocuments = db['form_options'].count();

// If no documents exist, insert default data into the collection
if (countDocuments === 0) {
  db['form_options'].insertMany([
    {
	"field": "Task",
	"name": "Clustering",
	"username": "default",
	"abbreviation": "CL"
},
{
	"field": "Task",
	"name": "Imputation",
	"username": "default",
	"abbreviation": "IM"
},
{
	"field": "Task",
	"name": "Marker gene identification",
	"username": "default",
	"abbreviation": "MGI"
},
{
	"field": "Task",
	"name": "Trajectory",
	"username": "default",
	"abbreviation": "TR"
},
{
	"field": "Task",
	"name": "Cell-cell communication",
	"username": "default",
	"abbreviation": "CCC"
},
{
	"field": "Task",
	"name": "Multi-omic data integration",
	"username": "default",
	"abbreviation": "MDI"
},
{
	"field": "Task",
	"name": "Gene regulatory relations",
	"username": "default",
	"abbreviation": "GRR"
},
{
	"field": "Task",
	"name": "Cell type identification",
	"username": "default",
	"abbreviation": "CTI"
},
{
	"field": "Task",
	"name": "Spatial",
	"username": "default",
	"abbreviation": "SP"
},
{
	"field": "Species",
	"name": "Human",
	"username": "default",
	"abbreviation": "h"
},
{
	"field": "Species",
	"name": "Mouse",
	"username": "default",
	"abbreviation": "m"
}
  ]);
  print("Default data inserted successfully.");
} else {
  print("Default data already exists in the collection. Skipping insertion.");
}

// Close the MongoDB connection
conn.close();
