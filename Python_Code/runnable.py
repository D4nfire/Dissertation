from intermine.webservice import Service
service = Service("http://yeastmine.yeastgenome.org/yeastmine/service", token = "YOUR-API-KEY")
query = service.new_query("Gene")
query.add_view(
    "primaryIdentifier", "secondaryIdentifier", "organism.shortName", "symbol",
    "name"
)
query.add_constraint("Gene", "IN", "systematic gene names", code = "A")

for row in query.rows():
    print(row["primaryIdentifier"], row["secondaryIdentifier"], row["organism.shortName"], \
        row["symbol"], row["name"])