-- MySQL dump 10.13  Distrib 5.1.41, for debian-linux-gnu (x86_64)
--
-- Host: db    Database: temp
-- ------------------------------------------------------
-- Server version	5.1.41-3ubuntu12.6-log

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `alignment`
--

DROP TABLE IF EXISTS `alignment`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `alignment` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `family_key` int(11) NOT NULL,
  `culled` tinyint(1) NOT NULL DEFAULT '0',
  `compressed` longblob NOT NULL,
  `format` enum('fasta','phylip','clustal','nexus') NOT NULL DEFAULT 'fasta',
  `filename` mediumtext,
  `timestamp` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  KEY `family_key` (`family_key`)
) ENGINE=MyISAM AUTO_INCREMENT=19139 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `alignment_colcull`
--

DROP TABLE IF EXISTS `alignment_colcull`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `alignment_colcull` (
  `alignment_key` int(11) NOT NULL,
  `column` smallint(5) unsigned NOT NULL,
  `gap_percentage` float unsigned NOT NULL,
  PRIMARY KEY (`alignment_key`,`column`),
  KEY `alignment_key` (`alignment_key`),
  KEY `column` (`column`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `alignment_protein`
--

DROP TABLE IF EXISTS `alignment_protein`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `alignment_protein` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `alignment_key` int(11) NOT NULL,
  `protein_key` int(11) DEFAULT NULL,
  `sequence` mediumtext,
  PRIMARY KEY (`id`),
  UNIQUE KEY `id` (`id`),
  KEY `sequence_key` (`protein_key`),
  KEY `alignment_key` (`alignment_key`),
  KEY `entry` (`protein_key`,`alignment_key`)
) ENGINE=MyISAM AUTO_INCREMENT=220195 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `alignment_seqcull`
--

DROP TABLE IF EXISTS `alignment_seqcull`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `alignment_seqcull` (
  `alignment_key` int(11) NOT NULL,
  `protein_key` int(11) unsigned NOT NULL,
  `score` float unsigned NOT NULL,
  PRIMARY KEY (`alignment_key`,`protein_key`),
  KEY `alignment_key` (`alignment_key`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `codeml`
--

DROP TABLE IF EXISTS `codeml`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `codeml` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `tree_key` int(11) NOT NULL,
  `model` int(11) NOT NULL,
  `filename` mediumtext,
  `timestamp` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  KEY `tree_key` (`tree_key`)
) ENGINE=MyISAM AUTO_INCREMENT=485 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `codeml_positive_selection`
--

DROP TABLE IF EXISTS `codeml_positive_selection`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `codeml_positive_selection` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `codeml_key` int(11) NOT NULL DEFAULT '0',
  `column` int(11) NOT NULL,
  `probability` double NOT NULL,
  `post_mean` double NOT NULL,
  `stderr` double DEFAULT NULL,
  PRIMARY KEY (`codeml_key`,`column`),
  UNIQUE KEY `id_2` (`id`),
  KEY `tree_key` (`codeml_key`),
  KEY `column` (`column`),
  KEY `probability` (`probability`),
  KEY `post_mean` (`post_mean`),
  KEY `stderr` (`stderr`),
  KEY `id` (`id`)
) ENGINE=MyISAM AUTO_INCREMENT=23949 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `family`
--

DROP TABLE IF EXISTS `family`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `family` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `oid_key` int(11) DEFAULT NULL,
  `name` varchar(150) NOT NULL,
  `description` tinytext,
  `manually_curated` tinyint(1) DEFAULT '0',
  `timestamp` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  UNIQUE KEY `name` (`name`),
  KEY `experiment_key` (`manually_curated`),
  FULLTEXT KEY `name_2` (`name`,`description`)
) ENGINE=MyISAM AUTO_INCREMENT=19715 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `family_protein`
--

DROP TABLE IF EXISTS `family_protein`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `family_protein` (
  `family_key` int(11) NOT NULL,
  `protein_key` int(11) NOT NULL DEFAULT '0',
  `seed` tinyint(1) DEFAULT '0',
  PRIMARY KEY (`family_key`,`protein_key`),
  KEY `sequence_key` (`protein_key`),
  KEY `seed` (`seed`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `oid`
--

DROP TABLE IF EXISTS `oid`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `oid` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `description` text,
  `timestamp` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=MyISAM AUTO_INCREMENT=4 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `oid_experiments`
--

DROP TABLE IF EXISTS `oid_experiments`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `oid_experiments` (
  `oid_key` int(11) NOT NULL,
  `experiment_key` int(11) NOT NULL,
  PRIMARY KEY (`oid_key`,`experiment_key`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tree`
--

DROP TABLE IF EXISTS `tree`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tree` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `alignment_key` int(11) NOT NULL,
  `compressed` longblob,
  `timestamp` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `filename` mediumtext,
  PRIMARY KEY (`id`),
  KEY `alignment_key` (`alignment_key`)
) ENGINE=MyISAM AUTO_INCREMENT=19136 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tree_node`
--

DROP TABLE IF EXISTS `tree_node`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tree_node` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `tree_key` int(11) NOT NULL,
  `protein_key` int(11) DEFAULT NULL,
  `ancestor_node` int(11) DEFAULT NULL,
  `branchlength` float DEFAULT '0',
  `branchlength_sum` float DEFAULT '0',
  PRIMARY KEY (`id`),
  KEY `id` (`id`),
  KEY `tree_key` (`tree_key`),
  KEY `sequence_key` (`protein_key`),
  KEY `entry` (`tree_key`,`protein_key`)
) ENGINE=MyISAM AUTO_INCREMENT=358134 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tree_node_diagchars`
--

DROP TABLE IF EXISTS `tree_node_diagchars`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tree_node_diagchars` (
  `tree_node_key` int(11) NOT NULL,
  `column` int(11) NOT NULL,
  `aa` mediumtext NOT NULL,
  PRIMARY KEY (`tree_node_key`,`column`),
  KEY `tree_node_key` (`tree_node_key`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2010-11-17 17:24:01
