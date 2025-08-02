<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<!--  Scicos.xsl  -->
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:strip-space elements="*"/>
	
<xsl:template match="/">
  <html>
	<head>
		<TITLE><xsl:value-of select="MAN/ModuleName"/></TITLE>
	</head>
	<body>
		<!--Table of Child-Links-->
		<A NAME="CHILD_LINKS">
			<STRONG>Subsections</STRONG>
		</A>
		<UL>
			<li/>
			<A HREF="#link1"><xsl:value-of select="MAN/links/link1"/></A>
			<li/>
			<A HREF="#link2"><xsl:value-of select="MAN/links/link2"/></A>
			<li/>
			<A HREF="#link3"><xsl:value-of select="MAN/links/link3"/></A>
			<li/>
			<A HREF="#link4"><xsl:value-of select="MAN/links/link4"/></A>
			<li/>
			<A HREF="#link5"><xsl:value-of select="MAN/links/link5"/></A>
			<li/>
			<A HREF="#link6"><xsl:value-of select="MAN/links/link6"/></A>
			<li/>
			<A HREF="#link5"><xsl:value-of select="MAN/links/link7"/></A>
			<li/>
			<A HREF="#link6"><xsl:value-of select="MAN/links/link8"/></A>
		</UL>
	<!--End of Table of Child-Links-->

		<hr/>
		<H2>
			<A NAME="link0"></A>
			<xsl:value-of select="MAN/ModuleName"/>
			<br/>
			<DIV  ALIGN="CENTER">
			  <xsl:element name="A">
				<xsl:element name="IMG">
				 <xsl:attribute name="src">           
				   <xsl:value-of select="MAN/ModuleIcon"/>
				 </xsl:attribute>
				</xsl:element>
			  </xsl:element>
			</DIV>
		</H2>
	
	<!--Internal Linkage-->
		<H3>
			<font color="blue">
				<A NAME="link1"><xsl:value-of select="MAN/links/link1"/></A>
			</font>
		</H3>
		<p><xsl:value-of select="MAN/Library"/></p>
		<br/>

		<H3>
			<font color="blue">
				<A NAME="link2"><xsl:value-of select="MAN/links/link2"/></A>
			</font>
		</H3>
		<PRE>
			<xsl:value-of select="MAN/Description"/>
		</PRE>
		<br/>

		<H3>
			<font color="blue">
				<A NAME="link3"><xsl:value-of select="MAN/links/link3"/></A>
			</font>
		</H3>
		<DIV ALIGN="CENTER">
			<TABLE BORDER="1">
			  <xsl:for-each select="MAN/DialogBoxes/Items/Item">
				<TR>
				  <TD ALIGN="CENTER"><xsl:value-of select="ItemTitle"/></TD>
				  <TD ALIGN="CENTER"><xsl:value-of select="ItemContent"/></TD>
				</TR>
			  </xsl:for-each>
			</TABLE>
		</DIV>
		<UL>
		  <xsl:for-each select="MAN/DialogBoxes/ItemExplain/Explain">
			  <LI>
				<xsl:value-of select="Subject"/>
				<xsl:value-of select="Content"/>
			  </LI>
		  </xsl:for-each>
		</UL>
		<br/>

		<H3>
			<font color="blue">
				<A NAME="link4"><xsl:value-of select="MAN/links/link4"/></A>
			</font>
		</H3>				
		<UL>
			<xsl:for-each select="MAN/DefaultProperties/Property">
				<LI>
					<xsl:value-of select="Subject"/>
					<b><xsl:value-of select="Content"/></b>
				</LI>
			</xsl:for-each>
		</UL>
		<br/>

		<H3>
			<font color="blue">
				<A NAME="link5"><xsl:value-of select="MAN/links/link5"/></A>
			</font>
		</H3>
		<H2>
			<DIV  ALIGN="CENTER">
			  <xsl:element name="A">
				<xsl:element name="IMG">
				  <xsl:attribute name="src">           
					<xsl:value-of select="MAN/ExampleDiagrams/ExampleIcon"/>
				  </xsl:attribute>
				</xsl:element>
			  </xsl:element>
			</DIV>
		</H2>
		<br/>

		<H3>
			<font color="blue">
				<A NAME="link6"><xsl:value-of select="MAN/links/link6"/></A>
			</font>
		</H3>
		<PRE>
			<xsl:value-of select="MAN/InterfacingFunction"/>
		</PRE>
		<br/>

		<H3>
			<font color="blue">
				<A NAME="link7"><xsl:value-of select="MAN/links/link7"/></A>
			</font>
		</H3>
		<PRE>
			<xsl:value-of select="MAN/ComputationalFunction"/>
		</PRE>
		<br/>

		<H3>
			<font color="blue">
				<A NAME="link8"><xsl:value-of select="MAN/links/link8"/></A>
			</font>
		</H3>
		  <p>
			  <xsl:for-each  select="MAN/SeeAlso/SEE_ALSO_ITEM">
				<xsl:element name="a">
					<xsl:attribute name="href">           
					<xsl:apply-templates select="./LINK/@PATH"/>
						<xsl:value-of select="./LINK"/>.htm</xsl:attribute>
					<xsl:apply-templates select="A|LINK|text()"/>
				</xsl:element>			
				<xsl:text disable-output-escaping="yes">,&amp;nbsp;</xsl:text>
				<xsl:text disable-output-escaping="yes">&amp;nbsp;</xsl:text>
			  </xsl:for-each>
		  </p>
		  <br/>
		
	<!--End of Internal Linkage-->
		
		<HR/>
		<author>
			<b><xsl:value-of select="MAN/authors"/></b>
			<i><xsl:value-of select="MAN/FilingDate"/></i>
		</author>
	</body>
  </html>
</xsl:template>
		
</xsl:stylesheet>
