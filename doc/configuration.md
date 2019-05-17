

## Database specification

The database is specified to commands using a URL as described below.

### MySQL

MySQL 
<pre>
mssql://[[USER][:PASSSPEC]@]HOST[:PORT]]/DBNAME
</pre>

The USER is the data user to use, HOST and PORT specified the port.
DBNAME is the name of the database.  The PASSSPEC element specified how
the password is obtained.

By default, the password is obtained on the local machine from
<pre>~/.my.cnf</pre> in the <pre>[client]</pre> section.
If PASSSPEC is a simple name, this is the section in <pre>~/.my.cnf</pre>
to access.  

The PASSSPEC may also be in the form <pre>MYCNF:SECTION</pre>.  Where MYCNF is
the local path to a <pre>MYCNF:SECTION</pre> configuration file. In this case,
<pre>/</pre> must be replaced with <pre>r%2F</pre>.  All of the strings parsed
from the are URL are URL decoded, so other encodings can be include.  Note that
if you use <pre>urllib.parse.quote()</pre> to generate the string, you need to
specify <pre>safe=''</pre> to escape forward-slashes.

### Sqlite3
<pre>
sqlite:///REL/PATH.db
sqlite:////ABS/PATH.db
</pre>

The <pre>sqlite</pre> URL specifies the locations of the file on a local file system.
The <pre>///</pre> form specifics a location relative to the current directory,
the <pre>////</pre> specifies an absolute path.
