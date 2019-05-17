"""
Unit tests for database support.
"""
import sys
import os
if __name__ == '__main__':
    rootDir = "../.."
    sys.path = [os.path.join(rootDir, "lib"),
                os.path.join(rootDir, "extern/pycbio/lib")] + sys.path
import unittest
from pycbio.sys.testCaseBase import TestCaseBase
from gencode_icedb.general.peeweeOps import DbUrl


class DbUrlTests(TestCaseBase):
    def testUrlParse1(self):
        u = DbUrl.parse("mysql://userA:pass@localhost:1234/db1")
        self.assertEqual(u,
                         DbUrl(url='mysql://userA:pass@localhost:1234/db1', scheme='mysql', netloc='userA:pass@localhost:1234', user='userA', passspec='pass', host='localhost', port='1234', database='db1'))

    def testUrlParse2(self):
        u = DbUrl.parse("mysql://userB@localhost:2345/db2")
        self.assertEqual(u,
                         DbUrl(url='mysql://userB@localhost:2345/db2', scheme='mysql', netloc='userB@localhost:2345', user='userB', passspec=None, host='localhost', port='2345', database='db2'))

    def testUrlParse3(self):
        u = DbUrl.parse("mysql://localhost:1234/db3")
        self.assertEqual(u,
                         DbUrl(url='mysql://localhost:1234/db3', scheme='mysql', netloc='localhost:1234', user=None, passspec=None, host='localhost', port='1234', database='db3'))

    def testUrlParse4(self):
        u = DbUrl.parse("mysql://localhost/db4")
        self.assertEqual(u,
                         DbUrl(url='mysql://localhost/db4', scheme='mysql', netloc='localhost', user=None, passspec=None, host='localhost', port=None, database='db4'))

    def testUrlParse5(self):
        u = DbUrl.parse("sqlite:///db5.db")
        self.assertEqual(u,
                         DbUrl(url='sqlite:///db5.db', scheme='sqlite', netloc=None, user=None, passspec=None, host=None, port=None, database='db5.db'))

    def testUrlParse6(self):
        u = DbUrl.parse("sqlite:////mnt/db6.db")
        self.assertEqual(u,
                         DbUrl(url='sqlite:////mnt/db6.db', scheme='sqlite', netloc=None, user=None, passspec=None, host=None, port=None, database='/mnt/db6.db'))

    def testUrlParse7(self):
        u = DbUrl.parse("mysql://dbuser:%7Efred%2F.my.cnf%3Agencode@localhost/gencodedb")
        self.assertEqual(u,
                         DbUrl(url='mysql://dbuser:%7Efred%2F.my.cnf%3Agencode@localhost/gencodedb', scheme='mysql', netloc='dbuser:~fred/.my.cnf:gencode@localhost', user='dbuser', passspec='~fred/.my.cnf:gencode', host='localhost', port=None, database='gencodedb'))

    def testUrlParse8(self):
        u = DbUrl.parse("mysql://dbuser:~fred%2F.my.cnf%3Agencode@localhost/gencodedb")
        self.assertEqual(u,
                         DbUrl(url='mysql://dbuser:~fred%2F.my.cnf%3Agencode@localhost/gencodedb', scheme='mysql', netloc='dbuser:~fred/.my.cnf:gencode@localhost', user='dbuser', passspec='~fred/.my.cnf:gencode', host='localhost', port=None, database='gencodedb'))


if __name__ == '__main__':
    unittest.main()
