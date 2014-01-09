
library(RMySQL);

mygocon = NULL;
hpfcon = NULL;
ddbcon = NULL;
hddbcon = NULL;

close_all_cons <- function()
{
	for(con in dbListConnections(MySQL())){mysqlCloseConnection(con)};
}

#kdrew: setup connection to mygo database
mygoSetup <- function()
{
	mygocon <<- dbConnect(MySQL(), group = "mygo");
	return(mygocon);
}

#kdrew: setup connection to hpf database
ddbSetup <- function()
{
	ddbcon <<- dbConnect(MySQL(), group = "ddb");
	return(ddbcon);
}

#kdrew: setup connection to hpf database
hpfSetup <- function()
{
	hpfcon <<- dbConnect(MySQL());#, group = "hpf");
	return(hpfcon);
}

#kdrew: setup connection to hpf database
hddbSetup <- function()
{
	hddbcon <<- dbConnect(MySQL(), group = "hddb");
	return(hddbcon);
}


#kdrew: qry the mygo database
mygo_qry <- function(statement)
{ 
	#kdrew: make connection if not already made
	if(is.null(mygocon))
	{ mygoSetup(); }
	qry(statement, mygocon); 
}

#kdrew: qry the ddb database
ddb_qry <- function(statement)
{ 
	#kdrew: make connection if not already made
	if(is.null(ddbcon))
	{ ddbSetup(); }
	qry(statement, ddbcon); 
}

#kdrew: qry the hpf database
hpf_qry <- function(statement)
{ 
	#kdrew: make connection if not already made
	if(is.null(hpfcon))
	{ hpfSetup(); }
	qry(statement, hpfcon); 
}

#kdrew: qry the hddb database
hddb_qry <- function(statement)
{ 
	#kdrew: make connection if not already made
	if(is.null(hpfcon))
	{ hddbSetup(); }
	qry(statement, hddbcon); 
}

#kdrew: shortcut for basic queries
qry <- function(statement, connection=con)
{
	dbGetQuery(connection, statement)
}

