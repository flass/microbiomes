
SET statement_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = off;
SET check_function_bodies = false;
SET client_min_messages = warning;
SET escape_string_warning = off;

--
-- Name: ltree; Type: EXTENSION; Owner: lassalle (needs SUPERUSER rights)
--
--~ 
-- unavailable in postgres 9.1 extensions ?
--~ SET search_path = public, pg_catalog;
--~ 
--~ CREATE EXTENSION ltree;

--
-- Name: taxonomy; Type: SCHEMA; Schema: -; Owner: lassalle
--

CREATE SCHEMA taxonomy;

ALTER SCHEMA taxonomy OWNER TO lassalle;



SET search_path = taxonomy, pg_catalog;


--
-- Name: citations; Type: TABLE; Schema: taxonomy; Owner: lassalle; Tablespace: 
--

CREATE TABLE citations (
        cit_id integer NOT NULL,
        cit_key VARCHAR(500),
        pubmed_id integer NOT NULL,
        medline_id integer NOT NULL,
        url text,
        text text,
        taxid_list text
);

ALTER TABLE taxonomy.citations OWNER TO lassalle;

--
-- Name: delnodes; Type: TABLE; Schema: taxonomy; Owner: lassalle; Tablespace: 
--

CREATE TABLE delnodes (
    tax_id integer NOT NULL
);


ALTER TABLE taxonomy.delnodes OWNER TO lassalle;

--
-- Name: division; Type: TABLE; Schema: taxonomy; Owner: lassalle; Tablespace: 
--

CREATE TABLE division (
    division_id integer NOT NULL,
    division_cde character(3),
    division_name text,
    comments text
);


ALTER TABLE taxonomy.division OWNER TO lassalle;

--
-- Name: gencode; Type: TABLE; Schema: taxonomy; Owner: lassalle; Tablespace: 
--

CREATE TABLE gencode (
    genetic_code_id character varying(3) NOT NULL,
    abbreviation text,
    name text,
    cde text,
    starts text
);


ALTER TABLE taxonomy.gencode OWNER TO lassalle;

--
-- Name: merged; Type: TABLE; Schema: taxonomy; Owner: lassalle; Tablespace: 
--

CREATE TABLE merged (
    old_tax_id integer NOT NULL,
    new_tax_id integer
);


ALTER TABLE taxonomy.merged OWNER TO lassalle;

--
-- Name: names; Type: TABLE; Schema: taxonomy; Owner: lassalle; Tablespace: 
--

CREATE TABLE names (
    tax_id integer,
    name_txt text,
    unique_name text,
    name_class text
);

CREATE INDEX ON taxonomy.names (tax_id);
CREATE INDEX ON taxonomy.names (name_class);
CREATE INDEX ON taxonomy.names (name_txt);

ALTER TABLE taxonomy.names OWNER TO lassalle;

--
-- Name: nodes; Type: TABLE; Schema: taxonomy; Owner: lassalle; Tablespace: 
--

-- stopped there


CREATE TABLE nodes (
    tax_id integer NOT NULL,
    parent_tax_id integer NOT NULL,
    rank character varying(50) NOT NULL,
    embl_code character(2),
    division_id integer,
    inherited_div_flag boolean,
    genetic_code_id character varying(3),
    inherited_gc_flag boolean,
    mitochondrial_genetic_code_id character varying(3),
    inherited_mgc_flag boolean,
    genbank_hidden_flag boolean,
    hidden_subtree_root_flag boolean,
    comments text
);


ALTER TABLE taxonomy.nodes OWNER TO lassalle;


CREATE TABLE gi_taxid_nucl (
	gi integer NOT NULL PRIMARY KEY,
    tax_id integer NOT NULL
);

ALTER TABLE taxonomy.gi_taxid_nucl OWNER TO lassalle;

CREATE INDEX ON taxonomy.gi_taxid_nucl (tax_id);


--
-- load data dumps;
--

\copy taxonomy.citations FROM '/home/florent/NCBI/Taxonomy/citations.dump'
\copy taxonomy.delnodes FROM '/home/florent/NCBI/Taxonomy/delnodes.dump'
\copy taxonomy.division FROM '/home/florent/NCBI/Taxonomy/division.dump'
\copy taxonomy.gencode FROM '/home/florent/NCBI/Taxonomy/gencode.dump'
\copy taxonomy.merged FROM '/home/florent/NCBI/Taxonomy/merged.dump'
\copy taxonomy.names FROM '/home/florent/NCBI/Taxonomy/names.dump'
\copy taxonomy.nodes FROM '/home/florent/NCBI/Taxonomy/nodes.dump'
\copy taxonomy.gi_taxid_nucl FROM '/home/florent/NCBI/Taxonomy/gi_taxid_nucl.dmp'

-- add path calumn to node table
--    path public.ltree NOT NULL


CREATE TABLE magda (
    name_txt text
);

ALTER TABLE taxonomy.magda OWNER TO lassalle;

\copy taxonomy.magda FROM '/home/florent/Desktop/colleagues/Magda/NewHostSpeciesListMagda2.modif'



CREATE TABLE phylosift_concat_reference (
    branch_id integer,
    tax_id integer
);

\copy taxonomy.phylosift_concat_reference FROM '/home/florent/oral_metagenomes/phylosift_v1.0.1/markers/concat.updated.annotated/concat.updated.taxonmap'
